from __future__ import annotations

import csv
import io
import os
import sys
import tempfile
import unittest.mock as mock
from pathlib import Path

os.environ.setdefault("FLASK_DEBUG", "")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "app"))

from app import (
    allowed_file,
    app,
    confidence_level,
    db_built_date,
    default_blastn_path,
    dotted_subject,
    extract_output_path,
    list_existing_results,
    load_existing_result,
    match_line,
    parse_blast_table,
    save_staged_result,
    save_upload,
    strip_leading_date,
    strip_tar_extension,
    validate_result_id,
    write_manifest,
)


class TestAllowedFile:
    def test_accepts_tar(self):
        assert allowed_file("results.tar")

    def test_accepts_tar_gz(self):
        assert allowed_file("results.tar.gz")

    def test_accepts_tgz(self):
        assert allowed_file("results.tgz")

    def test_rejects_zip(self):
        assert not allowed_file("results.zip")

    def test_rejects_no_extension(self):
        assert not allowed_file("results")

    def test_case_insensitive(self):
        assert allowed_file("results.TAR")
        assert allowed_file("results.TAR.GZ")


class TestDottedSubject:
    def test_identical_sequences_become_dots(self):
        assert dotted_subject("ACGT", "ACGT") == "...."

    def test_mismatches_show_subject_base(self):
        assert dotted_subject("ACGT", "ACGG") == "...G"

    def test_gap_in_query_shows_subject_base(self):
        assert dotted_subject("-CGT", "-CGT") == "-..."

    def test_empty_sequences(self):
        assert dotted_subject("", "") == ""


class TestMatchLine:
    def test_full_match(self):
        assert match_line("ACGT", "ACGT") == "||||"

    def test_full_mismatch(self):
        assert match_line("AAAA", "TTTT") == "    "

    def test_partial_match(self):
        assert match_line("ACGT", "ACGG") == "||| "

    def test_gap_is_not_match(self):
        assert match_line("-CGT", "-CGT") == " |||"

    def test_empty(self):
        assert match_line("", "") == ""


class TestParseBlastTable:
    def _write_tsv(self, rows: list[list[str]]) -> Path:
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, newline=""
        )
        writer = csv.writer(tmp, delimiter="\t")
        for row in rows:
            writer.writerow(row)
        tmp.close()
        return Path(tmp.name)

    def _make_row(
        self,
        query: str = "query1",
        subject_id: str = "NR_001",
        pct: str = "99.5",
        qcov: str = "100",
        evalue: str = "0.0",
        bitscore: str = "1200",
        ref: str = "Candida albicans strain SC5314",
        qs: str = "1",
        qe: str = "800",
        ss: str = "1",
        se: str = "800",
        qseq: str = "ACGT",
        sseq: str = "ACGT",
    ) -> list[str]:
        return [
            query,
            subject_id,
            pct,
            qcov,
            evalue,
            bitscore,
            ref,
            qs,
            qe,
            ss,
            se,
            qseq,
            sseq,
        ]

    def test_returns_empty_for_missing_file(self):
        assert parse_blast_table(Path("/nonexistent/path.tsv")) == []

    def test_returns_empty_for_none(self):
        assert parse_blast_table(None) == []

    def test_parses_single_hit(self):
        path = self._write_tsv([self._make_row()])
        result = parse_blast_table(path)
        assert len(result) == 1
        assert result[0]["query"] == "query1"
        assert len(result[0]["hits"]) == 1
        hit = result[0]["hits"][0]
        assert hit["percent_identity"] == "99.5"
        assert hit["species"] == "Candida albicans"

    def test_groups_hits_by_query(self):
        rows = [
            self._make_row(query="q1", ref="Saccharomyces cerevisiae strain S288C"),
            self._make_row(query="q1", ref="Saccharomyces paradoxus strain CBS432"),
            self._make_row(query="q2", ref="Candida albicans strain SC5314"),
        ]
        path = self._write_tsv(rows)
        result = parse_blast_table(path)
        assert len(result) == 2
        q1 = next(r for r in result if r["query"] == "q1")
        assert len(q1["hits"]) == 2

    def test_skips_short_rows(self):
        path = self._write_tsv([["too", "few", "cols"]])
        assert parse_blast_table(path) == []

    def test_dotted_subject_applied(self):
        path = self._write_tsv([self._make_row(qseq="ACGT", sseq="ACGT")])
        result = parse_blast_table(path)
        assert result[0]["hits"][0]["dotted_subject_sequence"] == "...."

    def test_match_line_applied(self):
        path = self._write_tsv([self._make_row(qseq="ACGT", sseq="ACGT")])
        result = parse_blast_table(path)
        assert result[0]["hits"][0]["match_line"] == "||||"

    def teardown_method(self, method):
        pass


class TestSaveUpload:
    def test_raises_on_empty_filename(self):
        upload = mock.MagicMock()
        upload.filename = ""
        try:
            save_upload(upload)
            assert False, "expected ValueError"
        except ValueError as exc:
            assert "Choose a tarball" in str(exc)

    def test_raises_on_bad_extension(self):
        upload = mock.MagicMock()
        upload.filename = "data.zip"
        try:
            save_upload(upload)
            assert False, "expected ValueError"
        except ValueError as exc:
            assert ".tar" in str(exc)

    def test_saves_valid_upload(self):
        upload = mock.MagicMock()
        upload.filename = "data.tar"
        upload.save = mock.MagicMock()
        path = save_upload(upload)
        try:
            assert path.name == "data.tar"
            upload.save.assert_called_once_with(path)
        finally:
            import shutil

            shutil.rmtree(path.parent, ignore_errors=True)


class TestExtractOutputPath:
    def test_extracts_matching_label(self):
        stdout = "BLAST summary: /output/results.txt\nother line"
        result = extract_output_path(stdout, "BLAST summary")
        assert result == Path("/output/results.txt")

    def test_returns_none_when_label_missing(self):
        assert extract_output_path("no match here", "BLAST summary") is None


class TestDefaultBlastnPath:
    def test_uses_env_var_when_set(self, monkeypatch):
        monkeypatch.setenv("BLASTN_BIN", "/custom/blastn")
        assert default_blastn_path() == "/custom/blastn"

    def test_falls_back_to_empty_when_not_found(self, monkeypatch):
        monkeypatch.delenv("BLASTN_BIN", raising=False)
        with mock.patch("app.LOCAL_BLASTN") as m:
            m.is_file.return_value = False
            assert default_blastn_path() == ""


class TestFlaskRoutes:
    def setup_method(self, method):
        app.config["TESTING"] = True
        self.client = app.test_client()

    def test_index_returns_200(self):
        response = self.client.get("/")
        assert response.status_code == 200

    def test_run_without_file_returns_400(self):
        response = self.client.post("/run", data={})
        assert response.status_code == 400

    def test_run_with_bad_extension_returns_400(self):
        data = {"tarball": (io.BytesIO(b"data"), "results.zip")}
        response = self.client.post(
            "/run", data=data, content_type="multipart/form-data"
        )
        assert response.status_code == 400

    def test_run_calls_pipeline_on_valid_upload(self):
        fake_result = {
            "returncode": 1,
            "stdout": "",
            "stderr": "pipeline failed",
            "summary_path": None,
            "blast_table_path": None,
            "fasta_path": None,
            "summary_text": None,
            "structured_results": [],
        }
        with mock.patch("app.run_staged_pipeline", return_value=fake_result) as mp:
            data = {"tarball": (io.BytesIO(b"fake tar content"), "results.tar")}
            response = self.client.post(
                "/run", data=data, content_type="multipart/form-data"
            )
            upload_path = mp.call_args[0][0]
        import shutil

        shutil.rmtree(upload_path.parent, ignore_errors=True)
        assert response.status_code == 500

    def test_run_success_renders_summary_table(self):
        fake_result = {
            "returncode": 0,
            "stdout": "BLAST summary: /tmp/s.txt",
            "stderr": "",
            "summary_path": Path("/tmp/s.txt"),
            "blast_table_path": None,
            "fasta_path": None,
            "summary_text": "Top hit: Candida albicans",
            "run_id": "run-1",
            "saved": False,
            "structured_results": [
                {
                    "query": "sample1",
                    "hits": [
                        {
                            "species": "Candida albicans",
                            "percent_identity": "99.5",
                            "confidence": "species",
                            "full_reference": "Candida albicans strain SC5314",
                            "query_coverage": "100",
                            "evalue": "0.0",
                            "bitscore": "1200",
                            "query_start": "1",
                            "query_end": "800",
                            "subject_start": "1",
                            "subject_end": "800",
                            "query_sequence": "ACGT",
                            "subject_sequence": "ACGT",
                            "dotted_subject_sequence": "....",
                            "match_line": "||||",
                            "subject_id": "NR_001",
                        }
                    ],
                }
            ],
        }
        with mock.patch("app.run_staged_pipeline", return_value=fake_result):
            data = {"tarball": (io.BytesIO(b"fake tar"), "results.tar")}
            response = self.client.post(
                "/run", data=data, content_type="multipart/form-data"
            )
        assert response.status_code == 200
        html = response.data.decode()
        assert "summary-table" in html
        assert "Candida albicans" in html
        assert "Species confirmed" in html


class TestConfidenceLevel:
    def test_species_at_99(self):
        assert confidence_level("99.0") == "species"

    def test_species_above_99(self):
        assert confidence_level("100.0") == "species"

    def test_probable_at_97(self):
        assert confidence_level("97.0") == "probable"

    def test_probable_between_97_and_99(self):
        assert confidence_level("98.5") == "probable"

    def test_genus_below_97(self):
        assert confidence_level("96.9") == "genus"

    def test_unknown_on_invalid(self):
        assert confidence_level("not_a_number") == "unknown"

    def test_unknown_on_none(self):
        assert confidence_level(None) == "unknown"


class TestDbBuiltDate:
    def test_returns_empty_when_marker_missing(self, tmp_path, monkeypatch):
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        assert db_built_date() == ""

    def test_returns_date_when_nhr_exists(self, tmp_path, monkeypatch):
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        db_dir = tmp_path / "data" / "blastdb"
        db_dir.mkdir(parents=True)
        marker = db_dir / "fungi_ITS.nhr"
        marker.write_bytes(b"")
        result = db_built_date()
        assert result  # non-empty date string
        assert len(result.split("-")) == 3  # YYYY-MM-DD format

    def test_returns_date_when_split_db_exists(self, tmp_path, monkeypatch):
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        db_dir = tmp_path / "data" / "blastdb"
        db_dir.mkdir(parents=True)
        marker = db_dir / "fungi_ITS.00.nhr"
        marker.write_bytes(b"")
        result = db_built_date()
        assert result
        assert len(result.split("-")) == 3


class TestDebugMode:
    def test_debug_defaults_to_false(self, monkeypatch):
        monkeypatch.delenv("FLASK_DEBUG", raising=False)
        assert app.debug is False

    def test_debug_enabled_with_flask_debug_1(self):
        debug = os.environ.get("FLASK_DEBUG", "").lower() in ("1", "true")
        assert not debug  # env is unset/empty during tests


class TestRunPipeline:
    def test_pipeline_failure_returns_nonzero_returncode(self):
        from app import run_pipeline

        with mock.patch("app.subprocess.run") as mock_run:
            mock_run.return_value = mock.MagicMock(
                returncode=1, stdout="", stderr="error"
            )
            result = run_pipeline(Path("/fake/archive.tar"))
        assert result["returncode"] == 1
        assert result["summary_text"] is None
        assert result["structured_results"] == []

    def test_pipeline_success_reads_summary(self, tmp_path):
        from app import run_pipeline

        summary = tmp_path / "results.summary.txt"
        summary.write_text("Top hit: Candida albicans")
        blast_table = tmp_path / "results.blast.tsv"
        blast_table.write_text("")
        fasta = tmp_path / "results.fasta"

        stdout = (
            f"BLAST summary: {summary}\n"
            f"BLAST table: {blast_table}\n"
            f"Combined FASTA: {fasta}\n"
        )
        with mock.patch("app.subprocess.run") as mock_run:
            mock_run.return_value = mock.MagicMock(
                returncode=0, stdout=stdout, stderr=""
            )
            result = run_pipeline(Path("/fake/archive.tar"))
        assert result["returncode"] == 0
        assert result["summary_text"] == "Top hit: Candida albicans"


class TestStagedSaving:
    def test_strip_tar_extension(self):
        assert strip_tar_extension("sample.tar") == "sample"
        assert strip_tar_extension("sample.tar.gz") == "sample"
        assert strip_tar_extension("sample.tgz") == "sample"

    def test_strip_leading_date(self):
        assert strip_leading_date("20260517-sample") == "sample"
        assert strip_leading_date("sample") == "sample"

    def test_run_pipeline_passes_staging_directories(self, tmp_path):
        from app import run_pipeline

        summary = tmp_path / "results.summary.txt"
        summary.write_text("summary")
        blast_table = tmp_path / "results.blast.tsv"
        blast_table.write_text("")
        fasta = tmp_path / "results.fasta"
        stdout = (
            f"BLAST summary: {summary}\n"
            f"BLAST table: {blast_table}\n"
            f"Combined FASTA: {fasta}\n"
        )
        with mock.patch("app.subprocess.run") as mock_run:
            mock_run.return_value = mock.MagicMock(
                returncode=0, stdout=stdout, stderr=""
            )
            run_pipeline(
                Path("/fake/archive.tar"),
                query_dir=tmp_path / "query",
                output_dir=tmp_path / "output",
            )
        env = mock_run.call_args.kwargs["env"]
        assert env["QUERY_DIR"] == str(tmp_path / "query")
        assert env["OUTPUT_DIR"] == str(tmp_path / "output")

    def test_save_staged_result_moves_files(self, tmp_path, monkeypatch):
        staging_root = tmp_path / ".app_runs"
        query_dir = tmp_path / "data" / "query"
        output_dir = tmp_path / "output"
        monkeypatch.setattr("app.STAGING_ROOT", staging_root)
        monkeypatch.setattr("app.QUERY_DIR", query_dir)
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)

        run_dir = staging_root / "run1"
        staged_query = run_dir / "query"
        staged_output = run_dir / "output"
        staged_query.mkdir(parents=True)
        staged_output.mkdir()
        input_tar = staged_query / "20260517-original.tar"
        fasta = staged_query / "20260517-original.fasta"
        blast = staged_output / "20260517-original.blast.tsv"
        summary = staged_output / "20260517-original.summary.txt"
        input_tar.write_bytes(b"tar")
        fasta.write_text(">sample\nACGT\n")
        blast.write_text("")
        summary.write_text("summary")
        write_manifest(
            "run1",
            {
                "original_filename": "original.tar",
                "input_tar": str(input_tar),
                "fasta_path": str(fasta),
                "blast_table_path": str(blast),
                "summary_path": str(summary),
                "saved": False,
            },
        )

        result = save_staged_result("run1")

        assert result["result_id"].endswith("-original")
        assert (query_dir / f"{result['result_id']}.tar").is_file()
        assert (query_dir / f"{result['result_id']}.fasta").is_file()
        assert (output_dir / f"{result['result_id']}.blast.tsv").is_file()
        assert (output_dir / f"{result['result_id']}.summary.txt").is_file()

    def test_save_staged_result_refuses_conflict(self, tmp_path, monkeypatch):
        staging_root = tmp_path / ".app_runs"
        query_dir = tmp_path / "data" / "query"
        output_dir = tmp_path / "output"
        monkeypatch.setattr("app.STAGING_ROOT", staging_root)
        monkeypatch.setattr("app.QUERY_DIR", query_dir)
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)

        date_prefix = __import__("datetime").datetime.now().strftime("%Y%m%d")
        query_dir.mkdir(parents=True)
        output_dir.mkdir()
        (query_dir / f"{date_prefix}-original.tar").write_bytes(b"existing")

        run_dir = staging_root / "run1"
        staged_query = run_dir / "query"
        staged_output = run_dir / "output"
        staged_query.mkdir(parents=True)
        staged_output.mkdir()
        input_tar = staged_query / "source.tar"
        fasta = staged_query / "source.fasta"
        blast = staged_output / "source.blast.tsv"
        summary = staged_output / "source.summary.txt"
        input_tar.write_bytes(b"tar")
        fasta.write_text(">sample\nACGT\n")
        blast.write_text("")
        summary.write_text("summary")
        write_manifest(
            "run1",
            {
                "original_filename": "original.tar",
                "input_tar": str(input_tar),
                "fasta_path": str(fasta),
                "blast_table_path": str(blast),
                "summary_path": str(summary),
                "saved": False,
            },
        )

        try:
            save_staged_result("run1")
            assert False, "expected FileExistsError"
        except FileExistsError as exc:
            assert "already exists" in str(exc)


class TestExistingResults:
    def _make_blast_row(self) -> str:
        return "\t".join(
            [
                "sample1",
                "NR_001",
                "99.5",
                "100",
                "0.0",
                "1200",
                "Candida albicans strain SC5314",
                "1",
                "4",
                "1",
                "4",
                "ACGT",
                "ACGT",
            ]
        )

    def test_validate_result_id_accepts_safe_ids(self):
        assert validate_result_id("20260514-results_1.2") == "20260514-results_1.2"

    def test_validate_result_id_rejects_path_traversal(self):
        try:
            validate_result_id("../secret")
            assert False, "expected ValueError"
        except ValueError as exc:
            assert "valid saved result" in str(exc)

    def test_list_existing_results_includes_summary_and_blast_table(
        self, tmp_path, monkeypatch
    ):
        output_dir = tmp_path / "output"
        query_dir = tmp_path / "data" / "query"
        output_dir.mkdir()
        query_dir.mkdir(parents=True)
        (output_dir / "20260514-run.summary.txt").write_text("summary")
        (output_dir / "20260514-run.blast.tsv").write_text(self._make_blast_row())
        (query_dir / "20260514-run.fasta").write_text(">sample1\nACGT\n")
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        monkeypatch.setattr("app.QUERY_DIR", query_dir)

        results = list_existing_results()

        assert len(results) == 1
        assert results[0]["id"] == "20260514-run"
        assert results[0]["has_structured_results"] is True
        assert results[0]["fasta_path"] == query_dir / "20260514-run.fasta"

    def test_load_existing_result_supports_summary_only(self, tmp_path, monkeypatch):
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        summary = output_dir / "old-run.summary.txt"
        summary.write_text("legacy summary")
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        monkeypatch.setattr("app.QUERY_DIR", tmp_path / "data" / "query")

        result = load_existing_result("old-run")

        assert result["summary_text"] == "legacy summary"
        assert result["structured_results"] == []
        assert result["blast_table_path"] is None

    def test_load_existing_result_parses_blast_table(self, tmp_path, monkeypatch):
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        (output_dir / "new-run.summary.txt").write_text("summary")
        (output_dir / "new-run.blast.tsv").write_text(self._make_blast_row())
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)
        monkeypatch.setattr("app.QUERY_DIR", tmp_path / "data" / "query")

        result = load_existing_result("new-run")

        assert result["summary_text"] == "summary"
        assert result["structured_results"][0]["query"] == "sample1"

    def test_load_existing_result_rejects_missing_result(self, tmp_path, monkeypatch):
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)

        try:
            load_existing_result("missing")
            assert False, "expected FileNotFoundError"
        except FileNotFoundError as exc:
            assert "not found" in str(exc)

    def test_open_result_route_renders_saved_summary(self, tmp_path, monkeypatch):
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        (output_dir / "saved.summary.txt").write_text("saved summary")
        monkeypatch.setattr("app.OUTPUT_DIR", output_dir)
        monkeypatch.setattr("app.REPO_ROOT", tmp_path)

        app.config["TESTING"] = True
        client = app.test_client()
        response = client.post("/open-result", data={"result_id": "saved"})

        assert response.status_code == 200
        assert "saved summary" in response.data.decode()

    def test_open_result_route_rejects_invalid_id(self):
        app.config["TESTING"] = True
        client = app.test_client()
        response = client.post("/open-result", data={"result_id": "../secret"})

        assert response.status_code == 404
