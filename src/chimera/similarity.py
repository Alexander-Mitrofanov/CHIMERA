"""Deterministic genome-similarity screening for leakage control.

The built-in engine uses canonical nucleotide k-mers and a bottom-k MinHash
sketch.  It is dependency-free and reproducible across Python versions because
it never uses Python's process-randomized ``hash`` function.  The reported
similarity is the Mash-style identity estimate derived from k-mer Jaccard; it
is an estimate, not an alignment ANI.  Publication workflows can instead pass
an explicit, versioned pairwise table produced by skani, FastANI, or another
validated method.
"""

from __future__ import annotations

import csv
import hashlib
import heapq
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from .errors import InputError
from .models import Genome, reverse_complement

SIMILARITY_TABLE_COLUMNS = (
    "query_genome_id",
    "reference_genome_id",
    "similarity",
    "coverage",
    "coverage_definition",
    "method",
)


def format_similarity_value(value: float | int | None) -> str:
    """Serialize fractional evidence without changing its binary-float value.

    Python's shortest-round-trip representation is deterministic and preserves
    values at scientific decision boundaries. Missing evidence is an empty TSV
    field, never a fabricated zero.
    """

    if value is None:
        return ""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError("similarity evidence must be a number or None")
    normalized = float(value)
    if not math.isfinite(normalized) or not 0.0 <= normalized <= 1.0:
        raise ValueError("similarity evidence must be finite and within [0, 1]")
    return repr(normalized)


def _canonical_kmer(kmer: str) -> str:
    reverse = reverse_complement(kmer)
    return kmer if kmer <= reverse else reverse


def _stable_hash64(value: str) -> int:
    digest = hashlib.blake2b(value.encode("ascii"), digest_size=8, person=b"CHIMERAk").digest()
    return int.from_bytes(digest, "big", signed=False)


@dataclass(frozen=True, slots=True)
class GenomeSketch:
    """A bottom-k sketch and enough state to identify exact sketches."""

    genome_id: str
    hashes: frozenset[int]
    k: int
    sketch_size: int
    unique_kmers: int
    genome_digest: str

    @property
    def exact(self) -> bool:
        """Whether every unique canonical k-mer is retained."""

        return self.unique_kmers <= self.sketch_size


@dataclass(frozen=True, slots=True)
class SimilarityHit:
    """Best observed training match for one query genome."""

    query_genome_id: str
    reference_genome_id: str | None
    similarity: float | None
    coverage: float | None
    method: str
    shared_hashes: int | None = None
    strict_gate_reference_genome_id: str | None = None
    strict_gate_similarity: float | None = None
    strict_gate_coverage: float | None = None
    strict_gate_method: str | None = None


def sketch_genome(genome: Genome, *, k: int = 21, sketch_size: int = 2_000) -> GenomeSketch:
    """Create a canonical, contig-boundary-aware genome sketch."""

    if k < 1:
        raise ValueError("k must be positive")
    if sketch_size < 1:
        raise ValueError("sketch_size must be positive")
    retained: set[int] = set()
    max_heap: list[int] = []
    overflow = False
    for contig in genome.contigs:
        sequence = contig.sequence
        for start in range(max(0, len(sequence) - k + 1)):
            kmer = sequence[start : start + k]
            # Ambiguity-bearing k-mers are uninformative for identity estimation.
            if any(base not in "ACGT" for base in kmer):
                continue
            hashed = _stable_hash64(_canonical_kmer(kmer))
            if hashed in retained:
                continue
            if len(retained) < sketch_size:
                retained.add(hashed)
                heapq.heappush(max_heap, -hashed)
                continue
            overflow = True
            current_largest = -max_heap[0]
            if hashed < current_largest:
                removed = -heapq.heapreplace(max_heap, -hashed)
                retained.remove(removed)
                retained.add(hashed)
    if not retained:
        raise InputError(
            f"Genome {genome.genome_id!r} has no unambiguous {k}-mers; "
            "choose a smaller --similarity-k or improve the reference sequence"
        )
    return GenomeSketch(
        genome_id=genome.genome_id,
        hashes=frozenset(retained),
        k=k,
        sketch_size=sketch_size,
        unique_kmers=sketch_size + 1 if overflow else len(retained),
        genome_digest=genome.digest,
    )


def sketch_jaccard(left: GenomeSketch, right: GenomeSketch) -> tuple[float, int]:
    """Estimate canonical k-mer Jaccard using a consistent bottom-k union sample."""

    if left.k != right.k:
        raise ValueError("Cannot compare sketches built with different k values")
    if left.genome_digest == right.genome_digest:
        return 1.0, min(len(left.hashes), len(right.hashes))
    if left.exact and right.exact:
        union = left.hashes | right.hashes
        shared = len(left.hashes & right.hashes)
        return (shared / len(union) if union else 0.0), shared

    sample_size = min(left.sketch_size, right.sketch_size)
    union_sample = sorted(left.hashes | right.hashes)[:sample_size]
    if not union_sample:
        return 0.0, 0
    shared = sum(value in left.hashes and value in right.hashes for value in union_sample)
    return shared / len(union_sample), shared


def mash_identity(jaccard: float, *, k: int) -> float | None:
    """Convert k-mer Jaccard to a Mash-style identity estimate.

    ``None`` denotes no shared canonical k-mer and is kept distinct from a
    measurable identity of zero.
    """

    if not 0.0 <= jaccard <= 1.0:
        raise ValueError("jaccard must be between 0 and 1")
    if jaccard == 0.0:
        return None
    if jaccard == 1.0:
        return 1.0
    distance = -(1.0 / k) * math.log((2.0 * jaccard) / (1.0 + jaccard))
    return max(0.0, min(1.0, 1.0 - distance))


def compare_sketches(query: GenomeSketch, reference: GenomeSketch) -> SimilarityHit:
    """Compare two sketches and return a self-describing result."""

    jaccard, shared = sketch_jaccard(query, reference)
    return SimilarityHit(
        query_genome_id=query.genome_id,
        reference_genome_id=reference.genome_id,
        similarity=mash_identity(jaccard, k=query.k),
        coverage=None,
        method=f"minhash-mash/k={query.k}/size={query.sketch_size}",
        shared_hashes=shared,
    )


def best_train_matches(
    train: Sequence[Genome],
    test: Sequence[Genome],
    *,
    k: int = 21,
    sketch_size: int = 2_000,
) -> dict[str, SimilarityHit]:
    """Compare every test genome with every training genome.

    Ties are resolved by the lexicographically smallest stable training ID, so
    record order cannot change the result.
    """

    if not train:
        raise ValueError("At least one training genome is required")
    train_sketches = [sketch_genome(genome, k=k, sketch_size=sketch_size) for genome in train]
    results: dict[str, SimilarityHit] = {}
    for genome in test:
        query = sketch_genome(genome, k=k, sketch_size=sketch_size)
        candidates = [compare_sketches(query, reference) for reference in train_sketches]
        best = min(
            candidates,
            key=lambda hit: (
                -(hit.similarity if hit.similarity is not None else -1.0),
                hit.reference_genome_id or "",
            ),
        )
        if best.similarity is None:
            best = replace(best, reference_genome_id=None)
        results[genome.genome_id] = best
    return results


def read_similarity_table(
    path: Path,
    *,
    query_ids: Iterable[str],
    reference_ids: Iterable[str],
    max_train_similarity: float = 0.95,
    min_similarity_coverage: float = 0.85,
) -> dict[str, SimilarityHit]:
    """Load an explicit pairwise table and select each query's best train hit.

    Detected similarity and aligned coverage must be fractions on [0, 1]. The
    table must contain exactly one row for every candidate-by-training pair,
    define coverage as ``aligned_fraction_shorter``, and identify a versioned
    method. A pair for which the external method detected no match is encoded
    by leaving *both* similarity and coverage blank; the pair row and its
    method provenance remain mandatory.

    The returned hit always describes the numerical maximum used for similarity
    stratification. Separate ``strict_gate_*`` fields preserve the strongest
    qualifying hit when it differs from that maximum.
    """

    for name, value in (
        ("max_train_similarity", max_train_similarity),
        ("min_similarity_coverage", min_similarity_coverage),
    ):
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be finite and within [0, 1]")

    query_sequence = tuple(query_ids)
    reference_sequence = tuple(reference_ids)
    expected_queries = set(query_sequence)
    expected_references = set(reference_sequence)
    if len(expected_queries) != len(query_sequence) or len(expected_references) != len(
        reference_sequence
    ):
        raise ValueError("query_ids and reference_ids must not contain duplicates")
    if not expected_queries or not expected_references:
        raise ValueError("query_ids and reference_ids must both be nonempty")
    by_query: dict[str, list[SimilarityHit]] = {query_id: [] for query_id in expected_queries}
    seen_pairs: set[tuple[str, str]] = set()
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as exc:
        raise InputError(f"Cannot read similarity table {path}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(SIMILARITY_TABLE_COLUMNS)
        missing_columns = required - set(reader.fieldnames or ())
        if missing_columns:
            raise InputError(
                f"Similarity table {path} is missing column(s): {', '.join(sorted(missing_columns))}"
            )
        for line_number, row in enumerate(reader, start=2):
            query_id = (row.get("query_genome_id") or "").strip()
            reference_id = (row.get("reference_genome_id") or "").strip() or None
            if query_id not in expected_queries:
                raise InputError(
                    f"{path}:{line_number}: unknown/non-test query_genome_id {query_id!r}"
                )
            if reference_id is None:
                raise InputError(
                    f"{path}:{line_number}: reference_genome_id is required in all-pairs mode"
                )
            if reference_id not in expected_references:
                raise InputError(
                    f"{path}:{line_number}: reference {reference_id!r} is not a training genome"
                )
            pair = (query_id, reference_id)
            if pair in seen_pairs:
                raise InputError(f"{path}:{line_number}: duplicate query/reference pair {pair!r}")
            seen_pairs.add(pair)
            raw_similarity = (row.get("similarity") or "").strip()
            raw_coverage = (row.get("coverage") or "").strip()
            if bool(raw_similarity) != bool(raw_coverage):
                raise InputError(
                    f"{path}:{line_number}: similarity and coverage must both be blank for "
                    "a nondetected pair, or both be numeric fractions on [0, 1]"
                )
            similarity_value: float | None = None
            coverage: float | None = None
            if raw_similarity:
                try:
                    similarity_value = float(raw_similarity)
                except ValueError as exc:
                    raise InputError(
                        f"{path}:{line_number}: similarity must be a number on [0, 1]"
                    ) from exc
                if not math.isfinite(similarity_value) or not 0.0 <= similarity_value <= 1.0:
                    raise InputError(
                        f"{path}:{line_number}: similarity {similarity_value!r} is outside [0, 1]"
                    )
                try:
                    coverage = float(raw_coverage)
                except ValueError as exc:
                    raise InputError(
                        f"{path}:{line_number}: coverage must be a number on [0, 1]"
                    ) from exc
                if not math.isfinite(coverage) or not 0.0 <= coverage <= 1.0:
                    raise InputError(
                        f"{path}:{line_number}: coverage {coverage!r} is outside [0, 1]"
                    )
            coverage_definition = (row.get("coverage_definition") or "").strip()
            if coverage_definition != "aligned_fraction_shorter":
                raise InputError(
                    f"{path}:{line_number}: coverage_definition must be 'aligned_fraction_shorter'"
                )
            method = (row.get("method") or "").strip()
            if not method or not any(character.isdigit() for character in method):
                raise InputError(
                    f"{path}:{line_number}: method must include a nonempty tool/version identifier"
                )
            by_query[query_id].append(
                SimilarityHit(
                    query_genome_id=query_id,
                    reference_genome_id=reference_id,
                    similarity=similarity_value,
                    coverage=coverage,
                    method=method,
                )
            )
    incomplete_queries = sorted(
        query for query, hits in by_query.items() if len(hits) != len(expected_references)
    )
    if incomplete_queries:
        preview = ", ".join(incomplete_queries[:5])
        suffix = " …" if len(incomplete_queries) > 5 else ""
        raise InputError(
            "Similarity table is not an exact candidate-by-training Cartesian matrix for "
            f"{len(incomplete_queries)} candidate(s): "
            f"{preview}{suffix}"
        )
    results: dict[str, SimilarityHit] = {}
    for query_id, hits in by_query.items():

        def crosses_strict_gate(hit: SimilarityHit) -> bool:
            return (
                hit.similarity is not None
                and hit.similarity > max_train_similarity
                and (hit.coverage is None or hit.coverage >= min_similarity_coverage)
            )

        detected_hits = [hit for hit in hits if hit.similarity is not None]
        if not detected_hits:
            provenance = min(hits, key=lambda hit: hit.reference_genome_id or "")
            results[query_id] = replace(provenance, reference_genome_id=None)
            continue
        maximum = min(
            detected_hits,
            key=lambda hit: (
                -(hit.similarity if hit.similarity is not None else -1.0),
                hit.reference_genome_id or "",
            ),
        )
        qualifying = [hit for hit in detected_hits if crosses_strict_gate(hit)]
        if qualifying:
            trigger = min(
                qualifying,
                key=lambda hit: (
                    -(hit.similarity if hit.similarity is not None else -1.0),
                    hit.reference_genome_id or "",
                ),
            )
            maximum = replace(
                maximum,
                strict_gate_reference_genome_id=trigger.reference_genome_id,
                strict_gate_similarity=trigger.similarity,
                strict_gate_coverage=trigger.coverage,
                strict_gate_method=trigger.method,
            )
        results[query_id] = maximum
    return results


def similarity_rows(hits: Mapping[str, SimilarityHit]) -> list[dict[str, object]]:
    """Return deterministic TSV-ready rows for provenance output."""

    return [
        {
            "query_genome_id": hit.query_genome_id,
            "reference_genome_id": hit.reference_genome_id or "",
            "similarity": format_similarity_value(hit.similarity),
            "coverage": format_similarity_value(hit.coverage),
            "method": hit.method,
        }
        for _, hit in sorted(hits.items())
    ]
