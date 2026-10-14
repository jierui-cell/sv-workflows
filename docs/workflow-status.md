# Method Scope

## Calling and Signatures

| Workflow | Output | Scope |
| --- | --- | --- |
| Manta, DELLY, SvABA | Paired somatic SV VCFs | WGS adapters; Manta also has capture mode |
| GRIDSS | Joint calls and native somatic subsets | Baseline without a panel of normals or BSgenome |
| BRASS | Annotated paired SVs | Validated GRCh37 resource bundle; missing statistics use labeled native defaults |
| COSMIC SV32 | Counts and SV1–SV12 exposures | GRCh37 WGS, caller-specific reference fitting |
| PCAWG/ClusterSV | On hold | No clustering/feature/NMF analysis; source and prior outputs retained |

Caller preflight checks declared GRCh37/GRCh38 reference dictionaries. BRASS's bundle adapter currently uses GRCh37. Signature routes use GRCh37 references.

## Assay Support

Select `sequencing_type` in caller YAML; use a separate run for each assay/capture design.

| Input | Implemented route |
| --- | --- |
| WGS | All five callers |
| WES | Manta with `--exome`; optional `target_bed` |
| Targeted sequencing | Manta with `--exome` and required `target_bed` |

[Manta's guide](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#exometargeted) explicitly supports exome/targeted mode and requires `--exome` even with call regions. The adapter sets it automatically.

This table describes our implemented adapters, not a claim that other tools cannot analyze capture reads. [SvABA](https://github.com/walaj/svaba/tree/v1.2.0) supports region-focused assembly, and [GRIDSS](https://github.com/PapenfussLab/gridss) discusses sparse targeted/exome data; equivalent assay-specific adapters remain future work. [DELLY's WGS-only restriction](https://github.com/dellytools/delly#cnv-calling) concerns its read-depth CNV command, which this SV workflow does not run.

Genome-wide signature methods use genome-wide breakpoint distributions. Capture designs change detection opportunity and the clustering null, so WES/panel callsets require a capture-aware analysis rather than pooling with WGS.

## Validation Coverage

The source workflows have processed paired WGS cohorts. The portable entry point has synthetic preflight/stub/resume tests; its real-pair equivalence run is not yet completed. Earlier ClusterSV/NMF numerical checks are retained as historical evidence, but those workflows are no longer executed in active CI. [Validation details](validation.md).

Slurm is compute-node tested. PBS/LSF/SGE are site templates; local execution works on workstations or allocated nodes. Public examples are synthetic and contain no patient data.

PCAWG/ClusterSV signature analysis and its CN-free extension will be revisited later. COSMIC SV32 and all five caller workflows remain active.
