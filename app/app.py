from __future__ import annotations

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


def default_blastn_path() -> str:
    if os.environ.get("BLASTN_BIN"):
        return os.environ["BLASTN_BIN"]
    if LOCAL_BLASTN.is_file():
        return str(LOCAL_BLASTN)
    return ""


def run_pipeline(upload_path: Path, blastn_path: str = "") -> dict[str, str | Path | None]:
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
    fasta_path = extract_output_path(completed.stdout, "Combined FASTA")
    summary_text = None

    if completed.returncode == 0:
        if summary_path is None or not summary_path.is_file():
            raise RuntimeError("The pipeline completed, but no BLAST summary file was found.")
        summary_text = summary_path.read_text(encoding="utf-8")

    return {
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "summary_path": summary_path,
        "fasta_path": fasta_path,
        "summary_text": summary_text,
    }


@app.get("/")
def index():
    return render_template("index.html", default_blastn_path=default_blastn_path())


@app.post("/run")
def run():
    upload = request.files.get("tarball")
    if upload is None:
        return render_template(
            "results.html",
            success=False,
            error="Choose a tarball before running the pipeline.",
        ), 400

    try:
        upload_path = save_upload(upload)
        blastn_path = request.form.get("blastn_path", "").strip()
        result = run_pipeline(upload_path, blastn_path)
    except Exception as exc:
        return render_template(
            "results.html",
            success=False,
            error=str(exc),
        ), 400

    success = result["returncode"] == 0
    error = None if success else "The pipeline did not finish successfully."

    return render_template(
        "results.html",
        success=success,
        error=error,
        stdout=result["stdout"],
        stderr=result["stderr"],
        summary_text=result["summary_text"],
        summary_path=result["summary_path"],
        fasta_path=result["fasta_path"],
        blastn_path=blastn_path,
        output_dir=OUTPUT_DIR,
    ), 200 if success else 500


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    app.run(host="127.0.0.1", port=port, debug=True)
