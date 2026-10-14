# Third-party software and resources

This repository contains workflow orchestration/adapters, not redistributed caller binaries, EGA data, model weights, reference bundles or licensed solvers. Each dependency retains its original license. Review upstream terms before use or redistribution.

| Dependency | Source / recorded identity |
| --- | --- |
| Nextflow | https://github.com/nextflow-io/nextflow ; tested 25.10.7 |
| Manta | https://github.com/Illumina/manta ; image digest in config/containers.tsv |
| DELLY | https://github.com/dellytools/delly ; image digest in config/containers.tsv |
| SvABA | https://github.com/walaj/svaba ; image digest in config/containers.tsv |
| GRIDSS | https://github.com/PapenfussLab/gridss ; image digest in config/containers.tsv |
| BRASS | https://github.com/cancerit/BRASS ; image digest in config/containers.tsv |
| Cue | https://github.com/PopicLab/cue ; user-supplied pinned code/model |
| LUMPY | https://github.com/arq5x/lumpy-sv ; discontinued historical tests only |
| JaBbA | https://github.com/mskilab-org/JaBbA ; licensed solver separate |
| ClusterSV | https://github.com/cancerit/ClusterSV ; AGPL-3.0-or-later, pinned commit in sv_signatures/pcawg/upstream.json |
| Panprostate feature code | https://github.com/panprostate/mutational-signatures ; pinned commit, fetched separately |
| SigProfiler / COSMIC | Python lockfile records tested packages; reference signatures retain provider terms |
| Replication timing | ENCODE/UW source lock in sv_signatures/pcawg/replication_sources.json |

The BRASS chromosome compatibility patch and adapters modify native sources only at runtime and record their hashes. Any redistribution of modified third-party source or combined work must comply with its upstream license, including applicable copyleft requirements. MIT licensing of our own independent workflow code does not replace these requirements.
