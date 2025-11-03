# jinja_template_nextflow

attempt to avoid boilerplate

testing capabilities of the free version...without giving it any code I care about

1. Generate nextflow pipeline

```
python generate_pipeline.py
```

2. Test it...

```
nextflow run generated_pipeline/main.nf \
  --bam_dir "data" \
  --reference "data/reference.fasta" \
  --callers "bcftools,freebayes"
```

Needs tuning...

```
 N E X T F L O W   ~  version 25.10.0

Launching `generated_pipeline/main.nf` [backstabbing_church] DSL2 - revision: ba5f24c6a1


Multi-Caller Variant Calling Pipeline
=====================================
BAM directory  : data
Reference      : data/reference.fasta
Output dir     : ./results
Callers        : bcftools,freebayes
Ensemble mode  : true
Min agreement  : 2

executor >  local (7)
[-        ] CALL_VARIANTS_GATK_HAPLOTYPECALLER -
[3d/5fc02d] CALL_VARIANTS_FREEBAYES (sample1)  [100%] 1 of 1 ✔
[15/f4467f] CALL_VARIANTS_BCFTOOLS (sample1)   [100%] 1 of 1 ✔
[-        ] CALL_VARIANTS_DEEPVARIANT          -
[-        ] CALL_VARIANTS_VARSCAN              -
[92/d159ae] NORMALIZE_VCF (sample1-freebayes)  [100%] 1 of 1 ✔
[e1/877187] COMPARE_CALLERS (sample1)          [  0%] 0 of 1
[ee/42f898] ENSEMBLE_CALLING (sample1)         [  0%] 0 of 1, retries: 2 ✘
[-        ] BENCHMARK_REPORT                   -

Pipeline completed!
===================
Status    : FAILED
Duration  : 9.8s
Results   : ./results

WARN: There's no process matching config selector: FILTER_SORT_BAM
[2a/1c1fdd] NOTE: Process `ENSEMBLE_CALLING (sample1)` terminated with an error exit status (1) -- Execution is retried (1)
[12/1d647a] NOTE: Process `ENSEMBLE_CALLING (sample1)` terminated with an error exit status (1) -- Execution is retried (2)
WARN: Killing running tasks (1)
ERROR ~ Error executing process > 'ENSEMBLE_CALLING (sample1)'

Caused by:
  Process `ENSEMBLE_CALLING (sample1)` terminated with an error exit status (1)

Command executed:

  # Merge all VCFs
  bcftools merge sample1.freebayes.norm.vcf.gz -Oz -o merged.vcf.gz
  tabix -p vcf merged.vcf.gz

  # Filter for variants called by at least N callers
  # This is simplified - production would use more sophisticated merging
  bcftools view merged.vcf.gz -Oz -o sample1.ensemble.vcf.gz

  echo "Ensemble calling complete" > sample1.ensemble_stats.txt
  echo "Minimum callers required: 2" >> sample1.ensemble_stats.txt
  echo "Total variants: $(bcftools view -H sample1.ensemble.vcf.gz | wc -l)" >> sample1.ensemble_stats.txt
  echo "Input callers: 1961" >> sample1.ensemble_stats.txt

Command exit status:
  1

Command output:
  (empty)

Command error:
  WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested

  About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.
           Note that only records from different files can be merged, never from the same file. For
           "vertical" merge take a look at "bcftools norm" instead.
  Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]

  Options:
          --force-samples               Resolve duplicate sample names
          --print-header                Print only the merged header and exit
          --use-header FILE             Use the provided header
      -0  --missing-to-ref              Assume genotypes at missing sites are 0/0
      -f, --apply-filters LIST          Require at least one of the listed FILTER strings (e.g. "PASS,.")
      -F, --filter-logic x|+            Remove filters if some input is PASS ("x"), or apply all filters ("+") [+]
      -g, --gvcf -|REF.FA               Merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max -M PL:max,AD:0
      -i, --info-rules TAG:METHOD,..    Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]
      -l, --file-list FILE              Read file names from the file
      -L, --local-alleles INT           EXPERIMENTAL: if more than <int> ALT alleles are encountered, drop FMT/PL and output LAA+LPL instead; 0=unlimited [0]
      -m, --merge STRING                Allow multiallelic records for <snps|indels|both|snp-ins-del|all|none|id>, see man page for details [both]
      -M, --missing-rules TAG:METHOD    Rules for replacing missing values in numeric vectors (.,0,max) when unknown allele <*> is not present [.]
          --no-index                    Merge unindexed files, the same chromosomal order is required and -r/-R are not allowed
          --no-version                  Do not append version and command line to the header
      -o, --output FILE                 Write output to a file [standard output]
      -O, --output-type u|b|v|z[0-9]    u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
      -r, --regions REGION              Restrict to comma-separated list of regions
      -R, --regions-file FILE           Restrict to regions listed in a file
          --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
          --threads INT                 Use multithreading with <int> worker threads [0]
          --write-index                 Automatically index the output files [off]

Work dir:
  /Users/jchang3/github/j23414/jinja_template_nextflow/work/ee/42f898e202704b9381a2844cb74f48

Container:
  quay.io/biocontainers/bcftools:1.18--h8b25389_0

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

3. Fix it in the main.nf, then fix it in the generate_pipeline.py

4. Think through the ramifications of this approach
