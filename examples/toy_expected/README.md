# Expected Synthetic Installation Output

After `python scripts/workflow.py test --profile test`, the default output is:

```text
results/toy/
├── PASS.txt           # Task execution/conversion assertions passed
├── toy.bedpe.tsv      # Three synthetic two-ended events
├── conversion.json   # Four VCF records collapsed to three retained events
└── audit.tsv          # Conversion decisions
```

The expected rows are in [toy.bedpe.tsv](toy.bedpe.tsv): a 20-kb deletion, a 30-kb tandem duplication and one interchromosomal adjacency. Reciprocal BND mates count once. These are synthetic fixtures, not patient data or a caller truth benchmark. Runtime JSON includes local paths and should not be committed.
