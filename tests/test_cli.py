from __future__ import annotations

import argparse
import json
import logging
import runpy
from pathlib import Path

import pytest

import chimera.cli as cli_module
from chimera.cli import _config_from_namespace, build_parser, legacy_main, main
from chimera.config import SplitKind
from chimera.errors import ConfigurationError, IntegrityError
from chimera.pipeline import BenchmarkResult


def write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    virus = tmp_path / "virus.fna"
    host = tmp_path / "host.fna"
    virus.write_text(">v1\n" + "ACGT" * 30 + "\n", encoding="utf-8")
    host.write_text(">h1\n" + "GATTACA" * 20 + "\n", encoding="utf-8")
    return virus, host


def test_root_help_exposes_recommended_one_click_command(capsys):
    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["--help"])
    assert raised.value.code == 0
    output = capsys.readouterr().out
    assert "suite" in output
    assert "Tests 2A-2E" in output
    assert "validate" in output


def test_suite_help_uses_plain_scientific_names(capsys):
    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["suite", "--help"])
    assert raised.value.code == 0
    output = capsys.readouterr().out
    assert "--release-date-cutoff" in output
    assert "--holdout-taxon" in output
    assert "--max-train-similarity" in output
    assert "--split" not in output
    assert "deposition" not in output.lower()

    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["generate", "--help"])
    assert raised.value.code == 0
    assert "--split" in capsys.readouterr().out


def test_cli_flags_override_config_and_resolve_from_cwd(tmp_path, monkeypatch):
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    config_path = config_dir / "suite.toml"
    config_path.write_text(
        """
[benchmark]
virus_paths = ["old-virus.fna"]
host_paths = ["old-host.fna"]
output_dir = "old-output"
test_fraction = 0.1
""".strip(),
        encoding="utf-8",
    )
    invocation = tmp_path / "invocation"
    invocation.mkdir()
    monkeypatch.chdir(invocation)
    args = build_parser().parse_args(
        [
            "generate",
            "--config",
            str(config_path),
            "--virus",
            "new-virus.fna",
            "--test-fraction",
            "0.3",
            "--split",
            "2b,2d",
        ]
    )

    config = _config_from_namespace(args)

    assert config.virus_paths == ((invocation / "new-virus.fna").resolve(),)
    assert config.host_paths == ((config_dir / "old-host.fna").resolve(),)
    assert config.test_fraction == 0.3
    assert config.splits == (SplitKind.GENOME, SplitKind.TEMPORAL)


def test_suite_rejects_a_selected_protocol_subset_from_toml(tmp_path):
    config_path = tmp_path / "subset.toml"
    config_path.write_text(
        '[benchmark]\nvirus_paths = ["v.fna"]\nhost_paths = ["h.fna"]\nsplits = ["genome"]\n',
        encoding="utf-8",
    )
    args = build_parser().parse_args(["suite", "--config", str(config_path)])

    with pytest.raises(ConfigurationError, match="suite always generates Tests 2A-2E"):
        _config_from_namespace(args)


def test_suite_parser_rejects_generate_only_split_flag(capsys):
    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["suite", "--split", "2b"])
    assert raised.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


def test_schema_command_is_machine_readable(capsys):
    assert main(["schema", "metadata"]) == 0
    header = capsys.readouterr().out.strip().split("\t")
    assert header[:2] == ["sequence_id", "genome_id"]
    assert "release_date" in header


def test_inspect_outputs_normalized_catalog(tmp_path, capsys):
    virus, host = write_inputs(tmp_path)
    result = main(
        [
            "inspect",
            "--virus",
            str(virus),
            "--host",
            str(host),
            "--json",
        ]
    )

    assert result == 0
    output = capsys.readouterr().out
    assert '"genome_id": "v1"' in output
    assert '"label": "host"' in output


def test_inspect_plain_text_contract(tmp_path, capsys):
    virus, host = write_inputs(tmp_path)

    assert main(["inspect", "--virus", str(virus), "--host", str(host)]) == 0

    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "genome_id\tlabel\tlength_nt\tcontigs\trelease_date\tsha256"
    assert len(lines) == 3


@pytest.mark.parametrize(
    ("flag", "value", "message"),
    [
        ("--test-fraction", "not-a-number", "must be a number"),
        ("--test-fraction", "-0.1", "must be a number"),
        ("--test-fraction", "1.1", "must be a number"),
        ("--fragment-length", "not-an-integer", "must be a positive integer"),
        ("--fragment-length", "0", "must be a positive integer"),
    ],
)
def test_numeric_flag_type_errors_are_actionable(flag, value, message, capsys):
    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["generate", flag, value])
    assert raised.value.code == 2
    assert message in capsys.readouterr().err


@pytest.mark.parametrize(
    "payload",
    [
        "[benchmark\n",
        "benchmark = []\n",
    ],
)
def test_invalid_toml_document_is_wrapped_as_configuration_error(tmp_path, payload):
    config_path = tmp_path / "invalid.toml"
    config_path.write_text(payload, encoding="utf-8")
    args = build_parser().parse_args(["generate", "--config", str(config_path)])

    with pytest.raises(ConfigurationError):
        _config_from_namespace(args)


def test_missing_toml_file_is_wrapped_as_configuration_error(tmp_path):
    args = build_parser().parse_args(["generate", "--config", str(tmp_path / "missing.toml")])
    with pytest.raises(ConfigurationError, match="Cannot read TOML"):
        _config_from_namespace(args)


def test_similarity_band_flags_merge_and_reject_non_table_config(tmp_path):
    config_path = tmp_path / "bands.toml"
    config_path.write_text(
        '[benchmark]\nvirus_paths = ["v.fna"]\nhost_paths = ["h.fna"]\n'
        "similarity_bands = { high = 0.9, moderate = 0.7, low = 0.3 }\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        ["generate", "--config", str(config_path), "--similarity-low", "0.2"]
    )
    assert _config_from_namespace(args).similarity_bands.low == 0.2

    config_path.write_text(
        '[benchmark]\nvirus_paths = ["v.fna"]\nhost_paths = ["h.fna"]\n'
        "similarity_bands = [0.9, 0.7, 0.3]\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        ["generate", "--config", str(config_path), "--similarity-low", "0.2"]
    )
    with pytest.raises(ConfigurationError, match="must be a table"):
        _config_from_namespace(args)


class _FakeValidationReport:
    def as_dict(self):
        return {"status": "pass", "checksums_verified": 3}

    def summary(self):
        return "validated fake bundle"


@pytest.mark.parametrize("json_output", [False, True])
def test_validate_cli_renders_human_and_json_reports(tmp_path, monkeypatch, capsys, json_output):
    monkeypatch.setattr("chimera.validation.validate_bundle", lambda _path: _FakeValidationReport())
    argv = ["validate", str(tmp_path)]
    if json_output:
        argv.append("--json")

    assert main(argv) == 0
    output = capsys.readouterr().out
    if json_output:
        assert json.loads(output)["checksums_verified"] == 3
    else:
        assert output.strip() == "PASS: validated fake bundle"


def test_validate_cli_maps_integrity_failure_to_exit_three(tmp_path, monkeypatch, caplog):
    def fail(_path):
        raise IntegrityError("tampered")

    monkeypatch.setattr("chimera.validation.validate_bundle", fail)
    with caplog.at_level(logging.ERROR, logger="chimera"):
        assert main(["validate", str(tmp_path)]) == 3
    assert "Integrity validation failed: tampered" in caplog.text


@pytest.mark.parametrize("dry_run", [False, True])
def test_generation_progress_branches_are_stable(tmp_path, monkeypatch, caplog, dry_run):
    result = BenchmarkResult(
        output_dir=tmp_path / "bundle",
        summary={"genomes": 4},
        dry_run=dry_run,
    )
    monkeypatch.setattr("chimera.pipeline.generate_benchmark", lambda *_args, **_kwargs: result)
    argv = [
        "generate",
        "--virus",
        str(tmp_path / "v.fna"),
        "--host",
        str(tmp_path / "h.fna"),
    ]
    if dry_run:
        argv.append("--dry-run")
    with caplog.at_level(logging.INFO, logger="chimera"):
        assert main(argv) == 0
    assert ("Preflight passed" if dry_run else "Publication bundle ready") in caplog.text


def test_main_maps_oserror_to_exit_two(monkeypatch, caplog):
    def fail(_args):
        raise OSError("unavailable")

    monkeypatch.setattr(cli_module, "_run_schema", fail)
    with caplog.at_level(logging.ERROR, logger="chimera"):
        assert main(["schema"]) == 2
    assert "unavailable" in caplog.text


def test_main_defensively_rejects_unknown_dispatched_command(monkeypatch):
    class FakeParser:
        def parse_args(self, _argv):
            return argparse.Namespace(command="unknown", quiet=False, log_level="INFO")

        def error(self, message):
            raise RuntimeError(message)

    monkeypatch.setattr(cli_module, "build_parser", lambda: FakeParser())
    with pytest.raises(RuntimeError, match="unsupported command"):
        main([])


def test_python_module_entry_point_delegates_to_cli(monkeypatch):
    monkeypatch.setattr(cli_module, "main", lambda: 23)
    imported = runpy.run_module("chimera.__main__", run_name="chimera.__main__.imported")

    with pytest.raises(SystemExit) as raised:
        runpy.run_module("chimera.__main__", run_name="__main__")

    assert imported["main"] is cli_module.main
    assert raised.value.code == 23


def test_legacy_entry_point_warns_and_delegates(monkeypatch, capsys):
    monkeypatch.setattr(cli_module, "main", lambda: 17)
    assert legacy_main() == 17
    assert "deprecated" in capsys.readouterr().err


def test_missing_required_sources_returns_documented_user_error(tmp_path):
    result = main(["suite", "--outdir", str(tmp_path / "out"), "--dry-run"])
    assert result == 2


def test_empty_outdir_flag_returns_exit_two_without_traceback(tmp_path, caplog, capsys):
    virus, host = write_inputs(tmp_path)

    with caplog.at_level(logging.ERROR, logger="chimera"):
        result = main(
            [
                "suite",
                "--virus",
                str(virus),
                "--host",
                str(host),
                "--outdir",
                "",
                "--dry-run",
            ]
        )

    captured = capsys.readouterr()
    assert result == 2
    assert "broad/protected directory" in caplog.text
    assert "Traceback" not in captured.err


@pytest.mark.parametrize(
    ("malformed_setting", "field_name"),
    [
        ("similarity_bands = [0.9, 0.7, 0.3]", "similarity_bands"),
        ("fragment_lengths = 31.5", "fragment_lengths"),
        ("holdout_taxa = 7", "holdout_taxa"),
        ('output_dir = ""', "output_dir"),
    ],
)
def test_malformed_toml_returns_exit_two_without_traceback(
    tmp_path, caplog, capsys, malformed_setting, field_name
):
    config_path = tmp_path / "malformed.toml"
    config_path.write_text(
        "\n".join(
            (
                "[benchmark]",
                'virus_paths = ["virus.fna"]',
                'host_paths = ["host.fna"]',
                malformed_setting,
            )
        ),
        encoding="utf-8",
    )

    with caplog.at_level(logging.ERROR, logger="chimera"):
        result = main(["suite", "--config", str(config_path), "--dry-run"])

    captured = capsys.readouterr()
    assert result == 2
    assert field_name in caplog.text
    assert "Traceback" not in captured.err
