# Tiny offline example

This deliberately synthetic fixture has two 320 nt viruses in two fictional
families and two 320 nt hosts. Each label has one genome released before and one
after the inclusive `2020-12-31` temporal cutoff. The sequences and taxonomy are
test data, not biological assertions.

From the repository root, after installing CHIMERA, run:

```console
chimera suite --config examples/tiny/chimera.toml
```

The command needs no network access and writes `examples/tiny/tiny-benchmark/`.
It creates all five protocols with eight fragments per genome, balanced between
31 nt and 61 nt. Re-running without `--force` refuses to replace the bundle;
add `--force` only when you intend to replace that recognized CHIMERA bundle.

Inspect inputs or check the completed bundle with:

```console
chimera inspect --virus examples/tiny/viruses.fna --host examples/tiny/hosts.fna --metadata examples/tiny/metadata.tsv
chimera validate examples/tiny/tiny-benchmark
```

For a no-write preflight, add `--dry-run` to the `suite` command. Paths inside
the TOML file are relative to that file, not to the shell's working directory.
