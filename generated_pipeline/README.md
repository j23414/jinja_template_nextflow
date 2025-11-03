# Multi-Caller Variant Calling Pipeline

Generated using Jinja2 templating.

## Included Callers

- **GATK HaplotypeCaller** (gatk_haplotypecaller)
- **FreeBayes** (freebayes)
- **BCFtools mpileup** (bcftools)
- **DeepVariant** (deepvariant)
- **VarScan** (varscan)

## Usage

```bash
# Run with all callers
nextflow run main.nf

# Run with specific callers
nextflow run main.nf --callers "gatk_haplotypecaller,freebayes"

# Disable ensemble calling
nextflow run main.nf --enable_ensemble false

# Custom parameters
nextflow run main.nf \
    --bam_dir /path/to/bams \
    --reference /path/to/ref.fasta \
    --outdir results \
    --min_callers_agree 3
```

## Output Structure

```
results/
├── variants/
│   ├── gatk_haplotypecaller/
│   ├── freebayes/
│   └── ...
├── comparison/
│   └── *_comparison.txt
├── ensemble/
│   └── *.ensemble.vcf.gz
└── reports/
    └── benchmark_report.html
```

## Adding New Callers

Edit `generate_pipeline.py` and add to `CALLERS_CONFIG`, then regenerate:

```python
python generate_pipeline.py
```
