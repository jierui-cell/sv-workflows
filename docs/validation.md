# Tests and Scientific Controls

Full-cohort launches include input/resource checks and per-pair preflight. COSMIC also runs its numerical clustering control automatically.

## Offline Tests

Inspect configured files without submitting jobs:

~~~bash
python scripts/workflow.py callers --config config/private-cohort.yaml --check-files
~~~

Run the synthetic test suite:

~~~bash
python scripts/test_installation.py
~~~

Tests cover configuration, sample roles, module setup, reference/region checks, VCF conversion, caller filters, automatic signature inputs and benchmark recording.

## Entry-Point Tests

On a workstation or allocated compute node:

~~~bash
python scripts/workflow.py test --profile test
python scripts/test_main_stub.py --out results/main_stub_check01
~~~

Both commands load the cohort YAML modules automatically. The toy test retains three SV events. The main-workflow test runs preflight on synthetic indexed BAMs, exercises four caller stubs and checks cached resume.

Add `--manta-image containers/images/manta-1.6.0--py27h9948957_6.sif` to include Manta, or `--with-exclusions` to test region handling. Use a fresh output directory for each test.

## Signature Controls

COSMIC’s clustering control checks clustered/nonclustered assignment and the pinned numerical behavior. PCAWG/ClusterSV analysis is paused.

## Optional Biological Evaluation

For accuracy studies, use normal-self controls, suitable positive truth sets, tumour/normal swaps, repeat runs and representative IGV review. Germline truth measures germline discovery; assess somatic accuracy with somatic truth or spike-ins.

These evaluations complement cohort analysis; they are not required launch steps.
