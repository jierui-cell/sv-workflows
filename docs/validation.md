# Tests and Scientific Controls

Setup checks are optional for users and run automatically in CI where dependencies are available. Full-cohort launches already include per-pair preflight and signature clustering self-tests.

## Offline Tests

To inspect your actual setup without submitting jobs:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
~~~

It checks nonempty/readable files for selected callers/pairs, including BAM indexes, reference indexes, images, BEDs and BRASS validation files. It does not stream BAMs, hash large inputs or execute callers. Real launches run the same check before Nextflow task submission. Unit tests cover missing files, empty images, broken links, caller/pair selection, standard filenames, relative paths, stable resume inputs and the no-submission failure path.

~~~bash
python scripts/test_installation.py
~~~

Synthetic tests cover YAML/schema validation, paired roles, option defaults, reference/index checks, BED regions, VCF orientation/mate collapse, filter adapters, missing BRASS statistics and benchmark recording.

Signature-input tests exercise automatic discovery for all five callers, native/compact outputs, caller/sample selection, completion records, hashed manifest creation, unchanged-input cache stability and the launcher handoff. File-check and preview modes submit no analysis tasks; manual external manifests remain supported.

## Entry-Point Tests

On a workstation or allocated node:

~~~bash
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

The toy task retains three SV events. The entry-point control runs actual `main.nf` preflight on indexed synthetic BAMs, exercises four caller stubs and verifies cached resume. Add `--manta-image containers/images/manta-1.6.0--py27h9948957_6.sif` for Manta configuration. `--with-exclusions` checks region policies. Choose a fresh test output directory.

## Signature Controls

PCAWG/ClusterSV numerical workflows are on hold and are not executed in current CI. Earlier controls and their results are retained for future review. Lightweight regression tests check that the paused route submits no tasks.

COSMIC's numerical clustering control checks actual clustered/nonclustered assignment. These controls test software and numerical behavior; they are separate from biological benchmarking.

## Optional Biological Evaluation

Use these when evaluating a new dataset, caller version or scientific configuration:

- Normal-self controls to examine retained somatic artifacts.
- Suitable positive truth sets and confident regions to measure precision/recall.
- Tumour/normal swaps and repeated-run comparisons to examine role handling and variation.
- Resolved-command review and representative read-level inspection in IGV.

GIAB germline truth measures germline discovery, not tumour-specific somatic precision. Keep benchmark scope explicit and use somatic truth/spike-ins for somatic accuracy.

These evaluations can be scheduled alongside cohort analysis; the launcher does not require a user-run pilot.
