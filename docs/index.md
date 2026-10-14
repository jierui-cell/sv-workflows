# SV Workflows

Paired tumour–normal SV calling with selectable signature analysis, designed for reproducible HPC research.

Start with [installation](setup.md), configure [your HPC](hpc.md), then validate [inputs and outputs](inputs-outputs.md). Scientific controls are explained in [validation](validation.md); signature choices and their limitations are in [signatures](signatures.md).

Five core callers are consolidated: Manta, DELLY, GRIDSS, SvABA and BRASS. Cue and JaBbA are research extensions; LUMPY is discontinued. COSMIC SV32 refitting and CN-free ClusterSV clustering are separate selectable workflows. No unavailable copy-number-dependent mechanism is fabricated.

This is a development release based on workflows previously used for paired tumour–normal WGS. The portable refactor still requires a real-pair equivalence pilot; successful workflow completion does not establish clinical accuracy. Public examples are synthetic; no patient data, cohort results or clinical metadata are distributed.

[Performance measurement](performance.md) · [Development and reproducibility](development.md) · [Source repository](https://github.com/jierui-cell/sv-workflows)
