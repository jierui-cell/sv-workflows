#!/usr/bin/env python3
"""Run the authors' fragCounter and CBS binaries/scripts; record exact commands."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from run_jabba import hash_file


def validate_options(settings):
    options = settings.get("coverage_options", {})
    if not isinstance(options, dict) or set(options) - {"window", "minmapq", "cnsignif"}:
        raise ValueError("Unsupported JaBbA coverage option")
    window = int(options.get("window", 1000))
    minmapq = int(options.get("minmapq", 20))
    cnsignif = float(options.get("cnsignif", 1e-5))
    if window < 100 or minmapq < 0 or minmapq > 60 or not 0 < cnsignif < 1:
        raise ValueError("Invalid fragCounter/CBS option")
    return window, minmapq, cnsignif


def bind_directories(project, bams):
    mounts = {project}
    mounts.update(bam.parent for bam in bams.values())
    return ",".join("{0}:{0}{1}".format(path, "" if path == project else ":ro")
                    for path in sorted(mounts, key=str))


def execute(command, log_path, commands_path, step):
    with commands_path.open("a") as record:
        record.write(json.dumps({"step": step, "argv": command}) + "\n")
    with log_path.open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("{} failed with exit {}; see {}".format(step, result.returncode, log_path))


def run(task):
    output = Path("output")
    output.mkdir(exist_ok=True)
    settings = task["settings"]
    window, minmapq, cnsignif = validate_options(settings)
    if settings.get("coverage_method") != "fragCounter-corrected-tumor-normal-ratio-1kb" or window != 1000:
        raise ValueError("Only the validated 1-kbp fragCounter coverage recipe is supported")
    if settings.get("cbs_source_commit") != "287c27648caa0ae05bb613e6adba7066f5ba17c2":
        raise ValueError("Unexpected CBS source revision")
    project = Path(task["project_dir"]).resolve()
    gcmap = Path(task["gcmap"]).resolve()
    frag_image = Path(task["frag_image"]).resolve()
    cbs_image = Path(task["cbs_image"]).resolve()
    cbs_script = Path(task["cbs_script"]).resolve()
    bams = {role: Path(task[role]).resolve() for role in ("tumor", "normal")}
    for path in (frag_image, cbs_image, cbs_script, *bams.values(),
                 gcmap / "gc1000.rds", gcmap / "map1000.rds"):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError("Missing/empty JaBbA coverage resource: " + str(path))
    commands = output / "commands.jsonl"
    commands.write_text("")
    # The BAMs need not live below the pipeline project directory. Mount their
    # actual parents read-only; the task work/output tree remains writable.
    binding = bind_directories(project, bams)
    for role in ("tumor", "normal"):
        temp = Path(role)
        temp.mkdir(exist_ok=True)
        command = ["singularity", "exec", "--cleanenv", "--bind", binding,
                   "--env", "R_DATATABLE_NUM_THREADS=1", str(frag_image),
                   "frag", "--bam", str(bams[role]), "--gcmapdir", str(gcmap),
                   "--window", str(window), "--minmapq", str(minmapq),
                   "--paired", "TRUE",
                   "--outdir", str(temp.resolve())]
        execute(command, output / (role + ".fragcounter.log"), commands, "fragcounter_" + role)
        cov = temp / "cov.rds"
        if not cov.is_file() or cov.stat().st_size == 0:
            raise RuntimeError("fragCounter did not produce " + str(cov))
        shutil.copy2(str(cov), str(output / (role + ".fragcounter.rds")))
    cbs_command = ["singularity", "exec", "--cleanenv", "--bind", binding,
                   "--env", "R_DATATABLE_NUM_THREADS=1", str(cbs_image),
                   "Rscript", str(cbs_script), "--tcov", str((Path("tumor") / "cov.rds").resolve()),
                   "--ncov", str((Path("normal") / "cov.rds").resolve()),
                   "--field", "reads.corrected", "--cnsignif", str(cnsignif),
                   "--name", task["pair_id"], "--outdir", str(output.resolve())]
    execute(cbs_command, output / "cbs.log", commands, "cbs_ratio")
    for filename in ("cov.rds", "seg.rds", "nseg.rds"):
        if not (output / filename).is_file() or (output / filename).stat().st_size == 0:
            raise RuntimeError("CBS ratio stage did not produce " + filename)
    (output / "cov.rds").rename(output / "coverage.rds")
    qc_script = ('x<-readRDS(commandArgs(TRUE)[1]); '
                 'stopifnot("ratio" %in% names(S4Vectors::mcols(x))); '
                 'z<-x$ratio; good<-is.finite(z)&z>0; '
                 'cat("bins\\tfinite_positive_bins\\tmedian_ratio\\n"); '
                 'cat(length(z),sum(good),median(z[good]),sep="\\t"); cat("\\n"); '
                 'stopifnot(sum(good)>0.5*length(z))')
    qc_command = ["singularity", "exec", "--cleanenv", "--bind", binding,
                  str(cbs_image), "Rscript", "-e", qc_script,
                  str((output / "coverage.rds").resolve())]
    execute(qc_command, output / "coverage.qc.tsv", commands, "coverage_qc")
    with (output / "artifacts.sha256").open("w") as handle:
        for path in (frag_image, cbs_image, cbs_script, gcmap / "gc1000.rds", gcmap / "map1000.rds"):
            handle.write("{}  {}\n".format(hash_file(path), path))
    (output / "effective-settings.json").write_text(json.dumps({
        "pair_id": task["pair_id"], "method": settings["coverage_method"],
        "settings": {"window": window, "minmapq": minmapq, "cnsignif": cnsignif},
        "cbs_source_commit": settings["cbs_source_commit"],
        "upstream": "https://github.com/mskilab-org/nf-gos"}, indent=2) + "\n")
    (output / "SUCCESS.json").write_text(json.dumps({"pair_id": task["pair_id"],
        "status": "success", "output": "coverage.rds"}, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: run_jabba_coverage.py task.json")
    run(json.loads(Path(sys.argv[1]).read_text()))
