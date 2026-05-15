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
    default_blastn_path,
    dotted_subject,
    extract_output_path,
    match_line,
    parse_blast_table,
    save_upload,
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
        assert match_line("ACGT", "ACGG") == "|||  "[:4]
        assert match_line("ACGT", "ACGG") == "|||" + " "

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
        assert path.name == "data.tar"
        upload.save.assert_called_once_with(path)


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
        with mock.patch("app.run_pipeline", return_value=fake_result):
            data = {"tarball": (io.BytesIO(b"fake tar content"), "results.tar")}
            response = self.client.post(
                "/run", data=data, content_type="multipart/form-data"
            )
        assert response.status_code == 500


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
