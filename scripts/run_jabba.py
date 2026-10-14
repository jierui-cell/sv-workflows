#!/usr/bin/env python3
"""Invoke upstream JaBbA in a pinned Singularity image; no SV algorithm is reimplemented."""

import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path


REQUIRED_OUTPUTS = ("jabba.simple.gg.rds", "jabba.simple.vcf", "jabba.seg", "opt.report.rds", "jabba.raw.rds")


def convergence_status(raw_gap, requested_gap, log_text):
    """Classify the native JaBbA result, not just its process exit code."""
    gap = float(raw_gap)
    target = float(requested_gap)
    if not math.isfinite(gap) or gap < 0:
        raise ValueError("JaBbA produced an invalid optimization gap")
    if not math.isfinite(target) or target <= 0:
        raise ValueError("Invalid requested JaBbA optimization gap")
    if "JaBbA did not converge" in log_text or gap > target:
        return {"status": "review_required_nonconverged", "converged": False,
                "final_epgap": gap, "requested_epgap": target}
    if "Done .. job output in:" not in log_text:
        raise ValueError("JaBbA log lacks its completion marker")
    return {"status": "converged", "converged": True,
            "final_epgap": gap, "requested_epgap": target}


def native_options(settings, purity, ploidy, cpus):
    solver = settings.get("solver", "gurobi")
    if solver not in ("cplex", "gurobi"):
        raise ValueError("JaBbA solver must be cplex or gurobi")
    options = settings.get("options", {})
    allowed = {"field", "ppmethod", "slack", "tilim", "epgap", "iterate", "cnsignif"}
    if not isinstance(options, dict) or set(options) - allowed:
        raise ValueError("Unknown JaBbA option(s): " + ", ".join(sorted(set(options) - allowed)))
    defaults = {"field": "ratio", "ppmethod": "ppgrid", "slack": 100,
                "tilim": 6000, "epgap": 1e-5, "iterate": 0, "cnsignif": 1e-5}
    defaults.update(options)
    if defaults["field"] not in ("ratio", "foreground") or defaults["ppmethod"] not in ("ppgrid", "ppurple", "sequenza"):
        raise ValueError("Invalid coverage field or purity/ploidy method")
    if (int(cpus) < 1 or float(settings["memory_gb"]) < 2 or
            float(defaults["slack"]) <= 0 or int(defaults["tilim"]) < 1 or
            float(defaults["epgap"]) <= 0 or int(defaults["iterate"]) < 0 or
            not 0 < float(defaults["cnsignif"]) < 1):
        raise ValueError("Invalid JaBbA compute or optimization setting")
    subgraph_cores = int(settings.get("subgraph_cores", 1))
    if subgraph_cores < 1 or subgraph_cores > int(cpus):
        raise ValueError("JaBbA subgraph_cores must be between 1 and allocated cpus")
    # This pinned upstream jba wrapper exposes --cores for subgraph forking but
    # leaves max.threads at its internal default of 12 per CPLEX solve.
    if "solver_threads" in settings:
        if int(settings["solver_threads"]) != 12:
            raise ValueError("Pinned JaBbA uses 12 solver threads; it does not expose a CLI override")
        if int(cpus) < 12 * subgraph_cores:
            raise ValueError("Allocate at least 12 CPUs per parallel JaBbA subgraph")
    args = ["--field", str(defaults["field"]), "--ppmethod", str(defaults["ppmethod"]),
            "--slack", str(defaults["slack"]), "--tilim", str(defaults["tilim"]),
            "--epgap", str(defaults["epgap"]), "--iterate", str(defaults["iterate"]),
            "--cnsignif", str(defaults["cnsignif"]), "--cores", str(subgraph_cores),
            "--mem", str(int(float(settings["memory_gb"]) * 0.8))]
    if solver == "gurobi":
        args.extend(["--gurobi", "TRUE"])
    if settings.get("use_metadata_purity_ploidy", False):
        if bool(purity) != bool(ploidy):
            raise ValueError("Purity and ploidy must both be supplied, or neither; no silent partial use")
        if purity and ploidy:
            args.extend(["--purity", str(purity), "--ploidy", str(ploidy)])
    return args, defaults


def hash_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(task):
    output = Path("output")
    output.mkdir(exist_ok=True)
    image = Path(task["image"]).resolve()
    junctions = Path(task["junctions"]).resolve()
    coverage = Path(task["coverage"]).resolve()
    for path in (image, junctions, coverage):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError("Missing or empty JaBbA input: " + str(path))
    settings = task["settings"]
    if settings.get("upstream_junction_caller") != "svaba":
        raise ValueError("Only the validated SvABA input policy is currently accepted")
    if not settings.get("coverage_method"):
        raise ValueError("Record the coverage preparation method in jabba.coverage_method")
    solver = settings.get("solver", "gurobi")
    if solver == "gurobi":
        home = os.environ.get("GUROBI_HOME")
        license_file = os.environ.get("GRB_LICENSE_FILE")
        if not home or not Path(home).is_dir() or not license_file or not Path(license_file).is_file():
            raise ValueError("Load a licensed Gurobi module in the Slurm task")
        solver_binds = {str(Path(home).resolve()), str(Path(license_file).resolve())}
        solver_env = ["--env", "GUROBI_HOME=" + home, "--env", "GRB_LICENSE_FILE=" + license_file,
                      "--env", "LD_LIBRARY_PATH=" + str(Path(home) / "lib")]
    elif solver == "cplex":
        cplex = Path(settings["cplex_dir"])
        if not cplex.is_absolute():
            cplex = Path(task["project_dir"]) / cplex
        cplex = cplex.resolve()
        if not (cplex / "cplex/include/ilcplex/cplex.h").is_file():
            raise ValueError("CPLEX_DIR lacks CPLEX headers: " + str(cplex))
        solver_binds = {str(cplex)}
        solver_env = ["--env", "CPLEX_DIR=" + str(cplex)]
    else:
        raise ValueError("JaBbA solver must be cplex or gurobi")
    r_lib = Path(settings["r_lib"])
    if not r_lib.is_absolute():
        r_lib = Path(task["project_dir"]) / r_lib
    r_lib = r_lib.resolve()
    if not r_lib.is_dir():
        raise ValueError("JaBbA Gurobi R library missing: " + str(r_lib))
    args, defaults = native_options(settings, task.get("purity"), task.get("ploidy"), settings["cpus"])
    filtered = output / "input.pass.bnd.vcf"
    filter_command = ["bcftools", "view", "-f", "PASS", "-i", 'INFO/SVTYPE="BND"',
                      "-Ov", "-o", str(filtered), str(junctions)]
    with (output / "commands.jsonl").open("w") as handle:
        handle.write(json.dumps({"step": "pass_bnd_filter", "argv": filter_command}) + "\n")
    subprocess.run(filter_command, check=True)
    records = 0
    ids = set()
    mates = []
    with filtered.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[6] != "PASS":
                raise ValueError("Malformed or non-PASS BND record after bcftools filtering")
            records += 1
            ids.add(fields[2])
            mates.extend(value[7:] for value in fields[7].split(";") if value.startswith("MATEID="))
    if set(mates) - ids:
        raise ValueError("Filtered SvABA BND VCF has orphan MATEID records")
    bind = ",".join(sorted({str(Path.cwd()), str(Path(task["project_dir"]).resolve()),
                            str(r_lib)} | solver_binds))
    prefix = ["singularity", "exec", "--cleanenv", "--bind", bind,
              "--env", "R_LIBS=" + str(r_lib),
              "--env", "R_LIBS_USER=" + str(r_lib)] + solver_env + [str(image)]
    locator = prefix + ["Rscript", "-e", "cat(system.file('extdata','jba',package='JaBbA'))"]
    location = subprocess.run(locator, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, check=True).stdout.strip()
    if not location:
        raise ValueError("Official JaBbA image does not contain an installed jba executable")
    command = prefix + ["Rscript", location, str(filtered.resolve()), str(coverage),
                        "--name", task["pair_id"], "--outdir", str(output.resolve())] + args
    (output / "effective-settings.json").write_text(json.dumps({
        "pair_id": task["pair_id"], "settings": settings, "native_defaults_and_overrides": defaults,
        "purity_ploidy_source": "samplesheet" if settings.get("use_metadata_purity_ploidy") and task.get("purity") and task.get("ploidy") else "JaBbA {} inference".format(defaults["ppmethod"]),
        "upstream_source": "https://github.com/mskilab-org/JaBbA",
        "junction_filter": "bcftools view -f PASS -i INFO/SVTYPE=BND",
        "pass_bnd_records": records,
        "solver": solver,
        "output_class": "graph-selected junctions dependent on SvABA input"}, indent=2) + "\n")
    with (output / "commands.jsonl").open("a") as handle:
        handle.write(json.dumps({"step": "jabba", "argv": command}) + "\n")
    with (output / "artifacts.sha256").open("w") as handle:
        for path in (image, junctions, coverage):
            handle.write("{}  {}\n".format(hash_file(path), path))
    with (output / "jabba.log").open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError("JaBbA exited {}; see output/jabba.log".format(completed.returncode))
    missing = [name for name in REQUIRED_OUTPUTS if not (output / name).is_file() or (output / name).stat().st_size == 0]
    if missing:
        raise RuntimeError("JaBbA exited zero but required output missing/empty: " + ", ".join(missing))
    with (output / "jabba.simple.vcf").open() as handle:
        if not handle.readline().startswith("##fileformat=VCF"):
            raise RuntimeError("JaBbA output VCF has no fileformat header")
    gap_command = prefix + ["Rscript", "-e", "x <- readRDS(commandArgs(TRUE)[1]); cat(x$epgap)",
                            str((output / "jabba.raw.rds").resolve())]
    gap = subprocess.check_output(gap_command, universal_newlines=True).strip()
    assessment = convergence_status(gap, defaults["epgap"], (output / "jabba.log").read_text())
    (output / "convergence.json").write_text(json.dumps(assessment, indent=2) + "\n")
    if not assessment["converged"]:
        (output / "REVIEW_REQUIRED.json").write_text(json.dumps({
            "pair_id": task["pair_id"], **assessment,
            "reason": "JaBbA returned a feasible graph but missed the configured optimization gap; do not treat as a validated callset"},
            indent=2) + "\n")
        return
    (output / "SUCCESS.json").write_text(json.dumps({"pair_id": task["pair_id"],
        "status": "success", "converged": True,
        "final_epgap": assessment["final_epgap"],
        "output_class": "graph-selected junctions, not independent de novo calls"}, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: run_jabba.py task.json")
    run(json.loads(Path(sys.argv[1]).read_text()))
