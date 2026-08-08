"""Deterministic fragment sampling and the Test 2A diagnostic split.

Sampling is performed independently for each genome and requested length.  A
coordinate is drawn uniformly from the union of all valid start coordinates on
eligible contigs, so longer contigs contribute in proportion to the number of
fragments they can produce.  Coordinates are sampled with replacement and can
therefore repeat, as they do in a shotgun library.

Random streams are derived from semantic BLAKE2 keys instead of from traversal
order.  Reordering genomes, contigs, or requested lengths consequently does not
change the generated dataset.  Truth coordinates are always zero-based and
half-open on the forward source contig, including for reverse-strand fragments.
"""

from __future__ import annotations

import hashlib
import math
import random
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from typing import Final, Literal, cast

from .errors import InputError, IntegrityError
from .models import Contig, Fragment, Genome, canonical_sequence_hash, reverse_complement

StrandMode = Literal["both", "forward"]

_SEED_PERSON: Final[bytes] = b"chimera.seed.v1"
_MAX_ATTEMPTS_PER_FRAGMENT: Final[int] = 10_000
_MIN_SPLIT_FRAGMENTS_PER_GENOME: Final[int] = 2
_REVERSE_STRAND_PROBABILITY: Final[float] = 0.5
_UNAMBIGUOUS_DNA: Final[frozenset[str]] = frozenset("ACGT")


def derive_seed(seed: int, *semantic_parts: str | int) -> int:
    """Derive a stable 128-bit sub-seed from an integer and semantic keys.

    Each component is type-tagged and length-prefixed, preventing ambiguous
    serializations such as ``("ab", "c")`` and ``("a", "bc")``.  The result
    does not depend on Python's process-randomized :func:`hash` function.

    Args:
        seed: Root seed for the benchmark run.  Negative integers are allowed.
        *semantic_parts: Stable names or integer indices describing the random
            stream, for example ``"coordinates"``, a genome ID, and a length.

    Returns:
        An unsigned integer suitable for constructing :class:`random.Random`.

    Raises:
        TypeError: If the seed or a semantic part has an unsupported type.
    """

    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer")

    digest = hashlib.blake2b(digest_size=16, person=_SEED_PERSON)
    for part in (seed, *semantic_parts):
        if isinstance(part, bool):
            raise TypeError("semantic seed parts must be strings or integers, not bool")
        if isinstance(part, int):
            encoded = str(part).encode("ascii")
            type_tag = b"i"
        elif isinstance(part, str):
            encoded = part.encode("utf-8")
            type_tag = b"s"
        else:
            raise TypeError(
                f"semantic seed parts must be strings or integers, got {type(part).__name__}"
            )
        digest.update(type_tag)
        digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
        digest.update(encoded)
    return int.from_bytes(digest.digest(), byteorder="big", signed=False)


def _validated_genomes(genomes: Iterable[Genome]) -> tuple[Genome, ...]:
    if isinstance(cast(object, genomes), (str, bytes)):
        raise InputError("genomes must be an iterable of Genome objects")
    try:
        materialized = tuple(genomes)
    except TypeError as error:
        raise InputError("genomes must be an iterable of Genome objects") from error
    if not materialized:
        raise InputError("at least one genome is required to generate fragments")
    if not all(isinstance(genome, Genome) for genome in materialized):
        raise InputError("genomes must contain only Genome objects")

    counts = Counter(genome.genome_id for genome in materialized)
    duplicates = sorted(genome_id for genome_id, count in counts.items() if count > 1)
    if duplicates:
        raise InputError(
            "genome IDs must be unique before fragment generation; duplicates: "
            + ", ".join(duplicates)
        )
    return tuple(sorted(materialized, key=lambda genome: (genome.genome_id, genome.digest)))


def _validated_lengths(fragment_lengths: Iterable[int]) -> tuple[int, ...]:
    if isinstance(fragment_lengths, (str, bytes)):
        raise InputError("fragment_lengths must be an iterable of positive integers")
    try:
        lengths = tuple(fragment_lengths)
    except TypeError as error:
        raise InputError("fragment_lengths must be an iterable of positive integers") from error
    if not lengths:
        raise InputError("fragment_lengths must contain at least one length")
    if any(isinstance(length, bool) or not isinstance(length, int) for length in lengths):
        raise InputError("fragment_lengths must contain only integers")
    invalid = sorted(length for length in lengths if length <= 0)
    if invalid:
        raise InputError(
            "fragment lengths must be positive; invalid value(s): "
            + ", ".join(str(length) for length in invalid)
        )
    duplicates = sorted(length for length, count in Counter(lengths).items() if count > 1)
    if duplicates:
        raise InputError(
            "fragment_lengths must not contain duplicates; duplicate value(s): "
            + ", ".join(str(length) for length in duplicates)
        )
    return tuple(sorted(lengths))


def _preflight_lengths(genomes: Sequence[Genome], lengths: Sequence[int]) -> None:
    unsupported: list[str] = []
    for genome in genomes:
        longest = max(contig.length for contig in genome.contigs)
        unsupported.extend(
            f"{genome.genome_id!r} cannot produce length {length} (longest contig: {longest})"
            for length in lengths
            if longest < length
        )
    if unsupported:
        details = "; ".join(unsupported)
        raise InputError(
            "Every genome must support every requested fragment length within one contig "
            f"(circular records may wrap their origin). Unsupported pair(s): {details}. "
            "Remove that length or provide a genome with a sufficiently long contig."
        )


def _balanced_length_counts(
    *,
    genome: Genome,
    lengths: Sequence[int],
    fragments_per_genome: int,
    seed: int,
) -> dict[int, int]:
    quotient, remainder = divmod(fragments_per_genome, len(lengths))
    counts = dict.fromkeys(lengths, quotient)
    remainder_order = list(lengths)
    random.Random(
        derive_seed(seed, "fragment-length-balance-v1", genome.genome_id, genome.digest)
    ).shuffle(remainder_order)
    for length in remainder_order[:remainder]:
        counts[length] += 1
    return counts


def _coordinate_space(
    genome: Genome,
    fragment_length: int,
) -> tuple[tuple[Contig, ...], tuple[int, ...]]:
    eligible = tuple(
        sorted(
            (contig for contig in genome.contigs if contig.length >= fragment_length),
            key=lambda contig: (contig.sequence_id, contig.digest),
        )
    )
    cumulative_ends: list[int] = []
    total = 0
    for contig in eligible:
        total += (
            contig.length if contig.topology == "circular" else contig.length - fragment_length + 1
        )
        cumulative_ends.append(total)
    return eligible, tuple(cumulative_ends)


def _draw_coordinate(
    contigs: Sequence[Contig],
    cumulative_ends: Sequence[int],
    rng: random.Random,
) -> tuple[Contig, int]:
    total_coordinates = cumulative_ends[-1]
    ticket = rng.randrange(total_coordinates)
    contig_index = bisect_right(cumulative_ends, ticket)
    previous_end = cumulative_ends[contig_index - 1] if contig_index else 0
    return contigs[contig_index], ticket - previous_end


def _ambiguous_fraction(sequence: str) -> float:
    return sum(base not in _UNAMBIGUOUS_DNA for base in sequence) / len(sequence)


def _opaque_fragment_id(
    *,
    seed: int,
    genome: Genome,
    contig: Contig,
    fragment_length: int,
    local_index: int,
    start: int,
    strand: str,
) -> str:
    identity = derive_seed(
        seed,
        "fragment-id-v1",
        genome.genome_id,
        genome.digest,
        contig.sequence_id,
        contig.digest,
        fragment_length,
        local_index,
        start,
        strand,
    )
    return f"frag-{identity:032x}"


def _validate_generation_options(
    *,
    fragments_per_genome: int,
    seed: int,
    strand_mode: str,
    max_ambiguous_fraction: float,
) -> float:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise InputError("seed must be an integer")
    if isinstance(fragments_per_genome, bool) or not isinstance(fragments_per_genome, int):
        raise InputError("fragments_per_genome must be an integer")
    if strand_mode not in ("both", "forward"):
        raise InputError("strand_mode must be 'both' or 'forward'")
    if isinstance(max_ambiguous_fraction, bool) or not isinstance(
        max_ambiguous_fraction, (int, float)
    ):
        raise InputError("max_ambiguous_fraction must be a finite number between 0 and 1")
    ambiguity_limit = float(max_ambiguous_fraction)
    if not math.isfinite(ambiguity_limit) or not 0.0 <= ambiguity_limit <= 1.0:
        raise InputError("max_ambiguous_fraction must be a finite number between 0 and 1")
    return ambiguity_limit


def generate_fragments(
    genomes: Iterable[Genome],
    *,
    fragment_lengths: Iterable[int],
    fragments_per_genome: int,
    seed: int,
    strand_mode: StrandMode = "both",
    max_ambiguous_fraction: float = 0.05,
) -> tuple[Fragment, ...]:
    """Generate balanced fixed-length fragments from labeled genomes.

    ``fragments_per_genome`` is the total for each genome, divided as evenly as
    possible among ``fragment_lengths``. It must provide at least two fragments
    per length, ensuring Test 2A can place every genome/length stratum on both
    sides. Start coordinates are sampled uniformly with
    replacement across every eligible coordinate on every contig.

    Args:
        genomes: Labeled source genomes.  Genome IDs must be unique.
        fragment_lengths: Distinct positive fixed lengths in nucleotides.
        fragments_per_genome: Exact total number emitted for each genome.
        seed: Root seed used to derive independent semantic random streams.
        strand_mode: ``"both"`` for equiprobable strands or ``"forward"``.
        max_ambiguous_fraction: Maximum fraction of non-ACGT IUPAC symbols in
            an emitted fragment, inclusive.

    Returns:
        An input-order-independent tuple of immutable :class:`Fragment`
        records with consecutive stable ordinals.

    Raises:
        InputError: If configuration or genomes are invalid, a genome cannot
            support every length, or bounded ambiguity rejection cannot fill
            the requested count.
    """

    ambiguity_limit = _validate_generation_options(
        fragments_per_genome=fragments_per_genome,
        seed=seed,
        strand_mode=strand_mode,
        max_ambiguous_fraction=max_ambiguous_fraction,
    )

    normalized_genomes = _validated_genomes(genomes)
    lengths = _validated_lengths(fragment_lengths)
    minimum_fragments = 2 * len(lengths)
    if fragments_per_genome < minimum_fragments:
        raise InputError(
            f"fragments_per_genome must be at least {minimum_fragments} so every genome "
            "contributes at least two fragments at every requested length"
        )
    _preflight_lengths(normalized_genomes, lengths)

    fragments: list[Fragment] = []
    seen_fragment_ids: set[str] = set()
    content_labels: dict[str, tuple[str, str]] = {}
    ordinal = 0
    for genome in normalized_genomes:
        counts = _balanced_length_counts(
            genome=genome,
            lengths=lengths,
            fragments_per_genome=fragments_per_genome,
            seed=seed,
        )
        for fragment_length in lengths:
            target_count = counts[fragment_length]
            contigs, cumulative_ends = _coordinate_space(genome, fragment_length)
            coordinate_rng = random.Random(
                derive_seed(
                    seed,
                    "fragment-coordinates-v1",
                    genome.genome_id,
                    genome.digest,
                    fragment_length,
                )
            )
            strand_rng = random.Random(
                derive_seed(
                    seed,
                    "fragment-strands-v1",
                    genome.genome_id,
                    genome.digest,
                    fragment_length,
                )
            )
            for local_index in range(target_count):
                for _attempt in range(_MAX_ATTEMPTS_PER_FRAGMENT):
                    contig, start = _draw_coordinate(contigs, cumulative_ends, coordinate_rng)
                    end = start + fragment_length
                    if end <= contig.length:
                        forward_sequence = contig.sequence[start:end]
                    else:
                        forward_sequence = (contig.sequence + contig.sequence)[start:end]
                    if _ambiguous_fraction(forward_sequence) <= ambiguity_limit:
                        break
                else:
                    total_coordinates = cumulative_ends[-1]
                    raise InputError(
                        f"Could not sample fragment {local_index + 1}/{target_count} of "
                        f"length {fragment_length} from genome {genome.genome_id!r} after "
                        f"{_MAX_ATTEMPTS_PER_FRAGMENT:,} random-coordinate attempts. No "
                        "sampled window satisfied "
                        f"max_ambiguous_fraction={ambiguity_limit:g} across "
                        f"{total_coordinates:,} eligible coordinate(s). Remove highly "
                        "ambiguous contigs, raise max_ambiguous_fraction, or request a "
                        "shorter fragment length."
                    )

                strand: Literal["+", "-"] = (
                    "-"
                    if strand_mode == "both" and strand_rng.random() < _REVERSE_STRAND_PROBABILITY
                    else "+"
                )
                sequence = (
                    reverse_complement(forward_sequence) if strand == "-" else forward_sequence
                )
                fragment_id = _opaque_fragment_id(
                    seed=seed,
                    genome=genome,
                    contig=contig,
                    fragment_length=fragment_length,
                    local_index=local_index,
                    start=start,
                    strand=strand,
                )
                if fragment_id in seen_fragment_ids:
                    raise IntegrityError(
                        "BLAKE2 fragment-ID collision detected; generation stopped before "
                        "writing an ambiguous truth record"
                    )
                seen_fragment_ids.add(fragment_id)
                content_digest = canonical_sequence_hash(sequence)
                previous = content_labels.get(content_digest)
                if previous is not None and previous[0] != genome.label.value:
                    raise InputError(
                        "Generated exact fragment content has contradictory virus/host labels: "
                        f"{previous[1]!r} ({previous[0]}) and {genome.genome_id!r} "
                        f"({genome.label.value}), canonical SHA-256 {content_digest}. "
                        "Curate integrated/proviral or contaminated sources before benchmarking."
                    )
                content_labels.setdefault(content_digest, (genome.label.value, genome.genome_id))
                fragments.append(
                    Fragment(
                        fragment_id=fragment_id,
                        sequence=sequence,
                        label=genome.label,
                        genome_id=genome.genome_id,
                        sequence_id=contig.sequence_id,
                        start=start,
                        end=end,
                        strand=strand,
                        ordinal=ordinal,
                    )
                )
                ordinal += 1
    return tuple(fragments)


def split_fragments_random(
    fragments: Iterable[Fragment],
    *,
    test_fraction: float,
    seed: int,
) -> tuple[tuple[Fragment, ...], tuple[Fragment, ...]]:
    """Create the Test 2A random-fragment diagnostic split.

    Membership is shuffled independently within every genome-by-length stratum
    so each source and requested length is represented in both partitions. The
    completed partitions are shuffled
    again with separate semantic streams, avoiding label- or source-grouped
    output order.  This diagnostic intentionally permits fragments from the
    same genome in train and test; genome-disjoint evaluation belongs to Test
    2B.

    Args:
        fragments: Generated fragments with globally unique fragment IDs.
        test_fraction: Desired per-genome/length test fraction, strictly between zero
            and one.  Counts are rounded to the nearest integer and clamped so
            each genome has at least one train and one test fragment.
        seed: Root seed used for membership and output-order streams.

    Returns:
        ``(train, test)`` tuples in deterministic shuffled order.

    Raises:
        InputError: If input is empty or invalid, IDs are duplicated, or a
            genome/length stratum has fewer than two fragments.
    """

    if isinstance(seed, bool) or not isinstance(seed, int):
        raise InputError("seed must be an integer")
    if isinstance(test_fraction, bool) or not isinstance(test_fraction, (int, float)):
        raise InputError("test_fraction must be a finite number between 0 and 1")
    normalized_fraction = float(test_fraction)
    if not math.isfinite(normalized_fraction) or not 0.0 < normalized_fraction < 1.0:
        raise InputError("test_fraction must be greater than 0 and less than 1")
    if isinstance(cast(object, fragments), (str, bytes)):
        raise InputError("fragments must be an iterable of Fragment objects")
    try:
        materialized = tuple(fragments)
    except TypeError as error:
        raise InputError("fragments must be an iterable of Fragment objects") from error
    if not materialized:
        raise InputError("at least one fragment is required for a random split")
    if not all(isinstance(fragment, Fragment) for fragment in materialized):
        raise InputError("fragments must contain only Fragment objects")

    id_counts = Counter(fragment.fragment_id for fragment in materialized)
    duplicate_ids = sorted(fragment_id for fragment_id, count in id_counts.items() if count > 1)
    if duplicate_ids:
        raise InputError(
            "fragment IDs must be unique before splitting; duplicates: " + ", ".join(duplicate_ids)
        )

    by_stratum: dict[tuple[str, int], list[Fragment]] = defaultdict(list)
    for fragment in materialized:
        by_stratum[(fragment.genome_id, fragment.length)].append(fragment)
    too_small = sorted(
        f"{genome_id}:{length}"
        for (genome_id, length), group in by_stratum.items()
        if len(group) < _MIN_SPLIT_FRAGMENTS_PER_GENOME
    )
    if too_small:
        raise InputError(
            "Test 2A requires at least two fragments in every genome/length stratum so "
            "both partitions are nonempty; insufficient strata: " + ", ".join(too_small)
        )

    train: list[Fragment] = []
    test: list[Fragment] = []
    for genome_id, length in sorted(by_stratum):
        group = sorted(by_stratum[(genome_id, length)], key=lambda fragment: fragment.fragment_id)
        random.Random(derive_seed(seed, "test-2a-membership-v2", genome_id, length)).shuffle(group)
        test_count = min(
            len(group) - 1,
            max(1, math.floor(len(group) * normalized_fraction + 0.5)),
        )
        test.extend(group[:test_count])
        train.extend(group[test_count:])

    train.sort(key=lambda fragment: fragment.fragment_id)
    test.sort(key=lambda fragment: fragment.fragment_id)
    random.Random(derive_seed(seed, "test-2a-train-order-v1")).shuffle(train)
    random.Random(derive_seed(seed, "test-2a-test-order-v1")).shuffle(test)

    train_ids = {fragment.fragment_id for fragment in train}
    test_ids = {fragment.fragment_id for fragment in test}
    if train_ids & test_ids:
        raise IntegrityError("Test 2A split produced overlapping fragment IDs")
    return tuple(train), tuple(test)


__all__ = [
    "StrandMode",
    "derive_seed",
    "generate_fragments",
    "split_fragments_random",
]
