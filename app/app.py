from __future__ import annotations

import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path

from flask import Flask, render_template, request
from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
PIPELINE_SCRIPT = REPO_ROOT / "script" / "02_process_seq_tarball_and_blast.sh"
OUTPUT_DIR = REPO_ROOT / "output"
ALLOWED_EXTENSIONS = (".tar", ".tar.gz", ".tgz")
LOCAL_BLASTN = Path("/usr/local/ncbi/blast/bin/blastn")
RESULT_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 128 * 1024 * 1024


def allowed_file(filename: str) -> bool:
    lower_name = filename.lower()
    return any(lower_name.endswith(extension) for extension in ALLOWED_EXTENSIONS)


def save_upload(upload: FileStorage) -> Path:
    filename = secure_filename(upload.filename or "")
    if not filename:
        raise ValueError("Choose a tarball before running the pipeline.")
    if not allowed_file(filename):
        raise ValueError("Upload a .tar, .tar.gz, or .tgz archive of .Seq files.")

    upload_dir = Path(tempfile.mkdtemp(prefix="its-upload-"))
    upload_path = upload_dir / filename
    upload.save(upload_path)
    return upload_path


def extract_output_path(stdout: str, label: str) -> Path | None:
    match = re.search(rf"^{re.escape(label)}:\s+(.+)$", stdout, re.MULTILINE)
    if not match:
        return None
    return Path(match.group(1).strip())


def db_built_date() -> str:
    db_dir = REPO_ROOT / "data" / "blastdb"
    marker = db_dir / "fungi_ITS.nhr"
    if not marker.is_file():
        marker = db_dir / "fungi_ITS.00.nhr"
    if not marker.is_file():
        return ""
    import datetime

    ts = marker.stat().st_mtime
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d")


# ITS species-level confidence thresholds (Irinyi et al. 2015, Med Mycol 53:313–37)
_SPECIES_THRESHOLD = 99.0
_PROBABLE_THRESHOLD = 97.0


def confidence_level(percent_identity: str) -> str:
    try:
        pct = float(percent_identity)
    except (ValueError, TypeError):
        return "unknown"
    if pct >= _SPECIES_THRESHOLD:
        return "species"
    if pct >= _PROBABLE_THRESHOLD:
        return "probable"
    return "genus"


def default_blastn_path() -> str:
    if os.environ.get("BLASTN_BIN"):
        return os.environ["BLASTN_BIN"]
    if LOCAL_BLASTN.is_file():
        return str(LOCAL_BLASTN)
    return ""


def dotted_subject(query_sequence: str, subject_sequence: str) -> str:
    dotted = []
    for query_base, subject_base in zip(query_sequence, subject_sequence):
        if query_base == subject_base and query_base != "-":
            dotted.append(".")
        else:
            dotted.append(subject_base)
    return "".join(dotted)


def match_line(query_sequence: str, subject_sequence: str) -> str:
    matches = []
    for query_base, subject_base in zip(query_sequence, subject_sequence):
        matches.append("|" if query_base == subject_base and query_base != "-" else " ")
    return "".join(matches)


def parse_blast_table(blast_table_path: Path | None) -> list[dict]:
    if blast_table_path is None or not blast_table_path.is_file():
        return []

    columns = (
        "query",
        "subject_id",
        "percent_identity",
        "query_coverage",
        "evalue",
        "bitscore",
        "full_reference",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "query_sequence",
        "subject_sequence",
    )

    queries: list[dict] = []
    by_query: dict[str, dict] = {}

    with blast_table_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < len(columns):
                continue

            hit = dict(zip(columns, row[: len(columns)]))
            species_words = hit["full_reference"].split()
            hit["species"] = " ".join(species_words[:2])
            hit["confidence"] = confidence_level(hit["percent_identity"])
            hit["dotted_subject_sequence"] = dotted_subject(
                hit["query_sequence"], hit["subject_sequence"]
            )
            hit["match_line"] = match_line(
                hit["query_sequence"], hit["subject_sequence"]
            )

            query = by_query.setdefault(
                hit["query"],
                {
                    "query": hit["query"],
                    "hits": [],
                },
            )
            if query not in queries:
                queries.append(query)
            query["hits"].append(hit)

    return queries


def list_existing_results() -> list[dict[str, object]]:
    results = []
    seen_ids = set()
    for summary_path in sorted(OUTPUT_DIR.glob("*.summary.txt"), reverse=True):
        result_id = summary_path.name.removesuffix(".summary.txt")
        seen_ids.add(result_id)
        blast_table_path = OUTPUT_DIR / f"{result_id}.blast.tsv"
        fasta_path = REPO_ROOT / "data" / "query" / f"{result_id}.fasta"
        results.append(
            {
                "id": result_id,
                "label": result_id,
                "summary_path": summary_path,
                "blast_table_path": blast_table_path if blast_table_path.is_file() else None,
                "fasta_path": fasta_path if fasta_path.is_file() else None,
                "has_structured_results": blast_table_path.is_file(),
            }
        )

    for blast_table_path in sorted(OUTPUT_DIR.glob("*.blast.tsv"), reverse=True):
        result_id = blast_table_path.name.removesuffix(".blast.tsv")
        if result_id in seen_ids:
            continue
        fasta_path = REPO_ROOT / "data" / "query" / f"{result_id}.fasta"
        results.append(
            {
                "id": result_id,
                "label": result_id,
                "summary_path": None,
                "blast_table_path": blast_table_path,
                "fasta_path": fasta_path if fasta_path.is_file() else None,
                "has_structured_results": True,
            }
        )

    return results


def validate_result_id(result_id: str) -> str:
    clean_id = result_id.strip()
    if not RESULT_ID_PATTERN.fullmatch(clean_id):
        raise ValueError("Choose a valid saved result.")
    return clean_id


def load_existing_result(result_id: str) -> dict[str, object]:
    clean_id = validate_result_id(result_id)
    summary_path = OUTPUT_DIR / f"{clean_id}.summary.txt"
    blast_table_path = OUTPUT_DIR / f"{clean_id}.blast.tsv"
    fasta_path = REPO_ROOT / "data" / "query" / f"{clean_id}.fasta"

    has_summary = summary_path.is_file()
    has_blast_table = blast_table_path.is_file()
    if not has_summary and not has_blast_table:
        raise FileNotFoundError("Saved result not found.")

    return {
        "returncode": 0,
        "stdout": "",
        "stderr": "",
        "summary_path": summary_path if has_summary else None,
        "blast_table_path": blast_table_path if has_blast_table else None,
        "fasta_path": fasta_path if fasta_path.is_file() else None,
        "summary_text": summary_path.read_text(encoding="utf-8") if has_summary else "",
        "structured_results": parse_blast_table(blast_table_path)
        if has_blast_table
        else [],
    }


def run_pipeline(upload_path: Path, blastn_path: str = "") -> dict[str, object]:
    env = os.environ.copy()
    if blastn_path:
        env["BLASTN_BIN"] = blastn_path

    completed = subprocess.run(
        ["bash", str(PIPELINE_SCRIPT), str(upload_path)],
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    summary_path = extract_output_path(completed.stdout, "BLAST summary")
    blast_table_path = extract_output_path(completed.stdout, "BLAST table")
    fasta_path = extract_output_path(completed.stdout, "Combined FASTA")
    summary_text = None
    structured_results = []

    if completed.returncode == 0:
        if summary_path is None or not summary_path.is_file():
            raise RuntimeError(
                "The pipeline completed, but no BLAST summary file was found."
            )
        summary_text = summary_path.read_text(encoding="utf-8")
        structured_results = parse_blast_table(blast_table_path)

    return {
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "summary_path": summary_path,
        "blast_table_path": blast_table_path,
        "fasta_path": fasta_path,
        "summary_text": summary_text,
        "structured_results": structured_results,
    }


@app.get("/")
def index():
    return render_template(
        "index.html",
        default_blastn_path=default_blastn_path(),
        db_built_date=db_built_date(),
        existing_results=list_existing_results(),
    )


@app.post("/run")
def run():
    upload = request.files.get("tarball")
    if upload is None:
        return (
            render_template(
                "results.html",
                success=False,
                error="Choose a tarball before running the pipeline.",
            ),
            400,
        )

    try:
        upload_path = save_upload(upload)
        blastn_path = request.form.get("blastn_path", "").strip()
        result = run_pipeline(upload_path, blastn_path)
    except Exception as exc:
        return (
            render_template(
                "results.html",
                success=False,
                error=str(exc),
            ),
            400,
        )

    success = result["returncode"] == 0
    error = None if success else "The pipeline did not finish successfully."

    return render_template(
        "results.html",
        success=success,
        error=error,
        stdout=result["stdout"],
        stderr=result["stderr"],
        summary_text=result["summary_text"],
        structured_results=result["structured_results"],
        summary_path=result["summary_path"],
        blast_table_path=result["blast_table_path"],
        fasta_path=result["fasta_path"],
        blastn_path=blastn_path,
        output_dir=OUTPUT_DIR,
        result_source="new",
    ), (200 if success else 500)


@app.post("/open-result")
def open_result():
    result_id = request.form.get("result_id", "")
    try:
        result = load_existing_result(result_id)
    except Exception as exc:
        return (
            render_template(
                "results.html",
                success=False,
                error=str(exc),
                stdout="",
                stderr="",
            ),
            404,
        )

    return render_template(
        "results.html",
        success=True,
        error=None,
        stdout="",
        stderr="",
        summary_text=result["summary_text"],
        structured_results=result["structured_results"],
        summary_path=result["summary_path"],
        blast_table_path=result["blast_table_path"],
        fasta_path=result["fasta_path"],
        blastn_path="",
        output_dir=OUTPUT_DIR,
        result_source="existing",
    )


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    debug = os.environ.get("FLASK_DEBUG", "").lower() in ("1", "true")
    app.run(host="127.0.0.1", port=port, debug=debug)
