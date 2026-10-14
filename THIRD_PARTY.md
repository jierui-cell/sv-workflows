# Third-party software and resources

This repository contains workflow orchestration/adapters, not redistributed caller binaries, EGA data, reference bundles or protected inputs. Each dependency retains its original license. Review upstream terms before use or redistribution.

| Dependency | Source / recorded identity |
| --- | --- |
| Nextflow | https://github.com/nextflow-io/nextflow ; tested 25.10.7 |
| Manta | https://github.com/Illumina/manta ; image digest in config/containers.tsv |
| DELLY | https://github.com/dellytools/delly ; image digest in config/containers.tsv |
| SvABA | https://github.com/walaj/svaba ; image digest in config/containers.tsv |
| GRIDSS | https://github.com/PapenfussLab/gridss ; image digest in config/containers.tsv |
| BRASS | https://github.com/cancerit/BRASS ; image digest in config/containers.tsv |
| ClusterSV | https://github.com/cancerit/ClusterSV ; AGPL-3.0-or-later, pinned commit in sv_signatures/pcawg/upstream.json |
| R NMF | https://github.com/renozao/NMF ; GPL (>=2), version 0.28 in sv_signatures/pcawg/renv.lock |
| SigProfiler / COSMIC | Python lockfile records tested packages; reference signatures retain provider terms |
| Replication timing | ENCODE/UW source lock in sv_signatures/pcawg/replication_sources.json |

The BRASS chromosome compatibility patch and adapters modify native sources only at runtime and record their hashes. Any redistribution of modified third-party source or combined work must comply with its upstream license, including applicable copyleft requirements. MIT licensing of our own independent workflow code does not replace these requirements.
