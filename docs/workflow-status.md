# Method Scope

## Calling and Signatures

| Workflow | Output and scope |
| --- | --- |
| Manta, DELLY, SvABA | Matched tumour–normal SV VCFs |
| GRIDSS | Joint calls and native somatic subsets; baseline without a panel of normals or BSgenome |
| BRASS | Annotated paired SVs using a GRCh37 resource bundle; missing purity/ploidy use recorded native defaults |
| COSMIC SV32 | GRCh37 WGS counts and SV1–SV12 reference-signature exposures |
| PCAWG/ClusterSV | On hold; code and prior outputs retained |

Caller preflight validates declared GRCh37/GRCh38 reference dictionaries. BRASS preparation and COSMIC SV32 currently use GRCh37.

## Assay Support

Set `sequencing_type` in YAML and use a separate run for each assay/capture design.

| Input | Implemented route |
| --- | --- |
| WGS | All five callers |
| WES | Manta `--exome`; optional `target_bed` |
| Targeted | Manta `--exome`; required `target_bed` |

The workflow sets [Manta’s capture mode](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#exometargeted) automatically. Other caller adapters currently use WGS settings. COSMIC SV32 uses genome-wide breakpoint distributions; capture callsets need a capture-aware signature model.

## Validation Coverage

Source workflows have processed paired WGS cohorts. CI tests the portable entry point with synthetic preflight, caller stubs, resume and COSMIC numerical controls; real-pair equivalence testing of the portable entry point remains to be completed. [Tests and controls](validation.md).

Slurm is compute-node tested. PBS/LSF/SGE profiles are configurable site templates. Public examples contain synthetic data only.
