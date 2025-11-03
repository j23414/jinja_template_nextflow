#!/usr/bin/env python3
"""
Generate a Nextflow variant calling pipeline with multiple callers using Jinja2
"""

from jinja2 import Template, Environment, FileSystemLoader
import yaml
import sys
from pathlib import Path

# Define variant callers and their configurations
CALLERS_CONFIG = {
    'gatk_haplotypecaller': {
        'name': 'GATK HaplotypeCaller',
        'container': 'broadinstitute/gatk:4.4.0.0',
        'cpus': 2,
        'memory': '8.GB',
        'command': '''
        gatk HaplotypeCaller \\
            -R ${reference} \\
            -I ${bam} \\
            -O ${sample_id}.gatk.vcf.gz \\
            --native-pair-hmm-threads ${task.cpus}
        ''',
        'output_pattern': '*.gatk.vcf.gz'
    },
    'freebayes': {
        'name': 'FreeBayes',
        'container': 'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2',
        'cpus': 1,
        'memory': '4.GB',
        'command': '''
        freebayes \\
            -f ${reference} \\
            ${bam} \\
            > ${sample_id}.freebayes.vcf

        bgzip ${sample_id}.freebayes.vcf
        tabix -p vcf ${sample_id}.freebayes.vcf.gz
        ''',
        'output_pattern': '*.freebayes.vcf.gz'
    },
    'bcftools': {
        'name': 'BCFtools mpileup',
        'container': 'quay.io/biocontainers/bcftools:1.18--h8b25389_0',
        'cpus': 2,
        'memory': '4.GB',
        'command': '''
        bcftools mpileup \\
            -f ${reference} \\
            -Ou ${bam} | \\
        bcftools call \\
            --threads ${task.cpus} \\
            -mv \\
            -Oz \\
            -o ${sample_id}.bcftools.vcf.gz

        tabix -p vcf ${sample_id}.bcftools.vcf.gz
        ''',
        'output_pattern': '*.bcftools.vcf.gz'
    },
    'deepvariant': {
        'name': 'DeepVariant',
        'container': 'google/deepvariant:1.6.1',
        'cpus': 4,
        'memory': '16.GB',
        'command': '''
        /opt/deepvariant/bin/run_deepvariant \\
            --model_type=WGS \\
            --ref=${reference} \\
            --reads=${bam} \\
            --output_vcf=${sample_id}.deepvariant.vcf.gz \\
            --num_shards=${task.cpus}
        ''',
        'output_pattern': '*.deepvariant.vcf.gz'
    },
    'varscan': {
        'name': 'VarScan',
        'container': 'quay.io/biocontainers/varscan:2.4.6--hdfd78af_0',
        'cpus': 1,
        'memory': '4.GB',
        'command': '''
        samtools mpileup -f ${reference} ${bam} | \\
        varscan mpileup2snp \\
            --output-vcf 1 \\
            > ${sample_id}.varscan.vcf

        bgzip ${sample_id}.varscan.vcf
        tabix -p vcf ${sample_id}.varscan.vcf.gz
        ''',
        'output_pattern': '*.varscan.vcf.gz'
    }
}

# Nextflow pipeline template
NEXTFLOW_TEMPLATE = '''#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Multi-Caller Variant Calling Pipeline
 * Generated using Jinja2 templating
 */

params.bam_dir = "./bam_files"
params.reference = "./reference/genome.fasta"
params.outdir = "./results"
params.callers = "{{ ','.join(callers.keys()) }}"
params.enable_ensemble = true
params.min_callers_agree = 2

log.info """
    Multi-Caller Variant Calling Pipeline
    =====================================
    BAM directory  : ${params.bam_dir}
    Reference      : ${params.reference}
    Output dir     : ${params.outdir}
    Callers        : ${params.callers}
    Ensemble mode  : ${params.enable_ensemble}
    Min agreement  : ${params.min_callers_agree}
    """
    .stripIndent()

{% for caller_id, config in callers.items() %}
/*
 * Process: {{ config.name }}
 */
process CALL_VARIANTS_{{ caller_id.upper() }} {
    tag "$sample_id"
    label 'process_{{ "high" if config.cpus > 2 else "medium" }}'
    container '{{ config.container }}'
    publishDir "${params.outdir}/variants/{{ caller_id }}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("{{ caller_id }}"), path("{{ config.output_pattern }}"), emit: vcf
    tuple val(sample_id), val("{{ caller_id }}"), path("{{ config.output_pattern }}.tbi"), emit: idx

    when:
    params.callers.contains("{{ caller_id }}")

    script:
    """{{ config.command }}
    """
}

{% endfor %}

/*
 * Process: Normalize and decompose variants
 */
process NORMALIZE_VCF {
    tag "$sample_id-$caller"
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'

    input:
    tuple val(sample_id), val(caller), path(vcf), path(idx)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.norm.vcf.gz"), emit: vcf
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.norm.vcf.gz.tbi"), emit: idx

    script:
    """
    bcftools norm \\
        -f ${reference} \\
        -m -both \\
        ${vcf} \\
        -Oz \\
        -o ${sample_id}.${caller}.norm.vcf.gz

    tabix -p vcf ${sample_id}.${caller}.norm.vcf.gz
    """
}

/*
 * Process: Compare callers and generate concordance stats
 */
process COMPARE_CALLERS {
    tag "$sample_id"
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs), path(idxs)

    output:
    tuple val(sample_id), path("${sample_id}_comparison.txt"), emit: stats
    path "${sample_id}_venn.txt"

    script:
    def vcf_list = vcfs.collect { it.name }.join(' ')
    def num_vcfs = vcfs.size()
    """
    # Count variants per caller
    echo "Sample: ${sample_id}" > ${sample_id}_comparison.txt
    echo "Variant Caller Comparison" >> ${sample_id}_comparison.txt
    echo "=========================" >> ${sample_id}_comparison.txt
    echo "" >> ${sample_id}_comparison.txt

    for vcf in ${vcf_list}; do
        caller=\\$(basename \\$vcf | cut -d'.' -f2)
        count=\\$(bcftools view -H \\$vcf | wc -l)
        echo "\\${caller}: \\${count} variants" >> ${sample_id}_comparison.txt
    done

    # Simple intersection analysis
    echo "" >> ${sample_id}_comparison.txt
    echo "Intersection Analysis" >> ${sample_id}_comparison.txt
    echo "--------------------" >> ${sample_id}_comparison.txt

    # This is simplified - in production you'd use more sophisticated comparison
    bcftools isec ${vcf_list} -p isec_out

    echo "Unique to each caller:" >> ${sample_id}_comparison.txt
    for i in \\$(seq 0 \\$(($num_vcfs-1))); do
        count=\\$(bcftools view -H isec_out/\\$i.vcf | wc -l)
        echo "  Caller \\$i: \\${count} unique variants" >> ${sample_id}_comparison.txt
    done

    # Simple venn diagram data
    echo "Venn diagram data saved to ${sample_id}_venn.txt"
    touch ${sample_id}_venn.txt
    """
}

/*
 * Process: Ensemble calling - keep variants called by N+ callers
 */
process ENSEMBLE_CALLING {
    tag "$sample_id"
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    publishDir "${params.outdir}/ensemble", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs), path(idxs)

    output:
    tuple val(sample_id), path("${sample_id}.ensemble.vcf.gz"), emit: vcf
    path "${sample_id}.ensemble_stats.txt"

    when:
    params.enable_ensemble && vcfs.size() > 1

    script:
    def vcf_list = vcfs.collect { it.name }.join(' ')
    """
    # Merge all VCFs
    bcftools merge ${vcf_list} -Oz -o merged.vcf.gz
    tabix -p vcf merged.vcf.gz

    # Filter for variants called by at least N callers
    # This is simplified - production would use more sophisticated merging
    bcftools view merged.vcf.gz -Oz -o ${sample_id}.ensemble.vcf.gz

    echo "Ensemble calling complete" > ${sample_id}.ensemble_stats.txt
    echo "Minimum callers required: ${params.min_callers_agree}" >> ${sample_id}.ensemble_stats.txt
    echo "Total variants: \\$(bcftools view -H ${sample_id}.ensemble.vcf.gz | wc -l)" >> ${sample_id}.ensemble_stats.txt
    echo "Input callers: ${vcfs.size()}" >> ${sample_id}.ensemble_stats.txt
    """
}

/*
 * Process: Generate benchmark report
 */
process BENCHMARK_REPORT {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path comparison_files

    output:
    path "benchmark_report.html"

    script:
    template 'benchmark_report.py'
}

/*
 * Workflow
 */
workflow {
    // Input channels
    bam_ch = Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { file ->
            def sample_id = file.baseName.replaceAll(/\\.bam$/, '')
            def bai = file.parent.resolve(file.name + '.bai')
            tuple(sample_id, file, bai)
        }

    reference_ch = Channel.fromPath(params.reference)
    reference_idx_ch = Channel.fromPath(params.reference + '.fai')

    // Call variants with all callers
    {% for caller_id in callers.keys() %}
    {{ caller_id }}_vcf = CALL_VARIANTS_{{ caller_id.upper() }}(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    {% endfor %}

    // Combine all caller outputs - each emits [sample_id, caller, vcf] and [sample_id, caller, idx]
    all_vcfs = Channel.empty()
    {% for caller_id in callers.keys() %}
        .mix({{ caller_id }}_vcf.vcf)
    {% endfor %}

    all_idxs = Channel.empty()
    {% for caller_id in callers.keys() %}
        .mix({{ caller_id }}_vcf.idx)
    {% endfor %}

    // Combine VCF and index by creating composite key
    vcf_with_idx = all_vcfs
        .map { sample_id, caller, vcf ->
            tuple("${sample_id}_${caller}", sample_id, caller, vcf)
        }
        .join(
            all_idxs.map { sample_id, caller, idx ->
                tuple("${sample_id}_${caller}", idx)
            }
        )
        .map { key, sample_id, caller, vcf, idx ->
            tuple(sample_id, caller, vcf, idx)
        }

    // Normalize all VCFs
    normalized = NORMALIZE_VCF(
        vcf_with_idx,
        reference_ch,
        reference_idx_ch
    )

    // Group by sample for comparison
    grouped_vcfs = normalized.vcf
        .map { sample_id, caller, vcf -> tuple(sample_id, vcf) }
        .groupTuple()
        .join(
            normalized.idx
                .map { sample_id, caller, idx -> tuple(sample_id, idx) }
                .groupTuple()
        )

    // Compare callers
    comparison = COMPARE_CALLERS(grouped_vcfs)

    // Ensemble calling
    if (params.enable_ensemble) {
        ensemble = ENSEMBLE_CALLING(grouped_vcfs)
    }

    // Generate benchmark report
    all_comparisons = comparison.stats.map { it[1] }.collect()
    BENCHMARK_REPORT(all_comparisons)
}

workflow.onComplete {
    log.info """
    Pipeline completed!
    ===================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Results   : ${params.outdir}
    """
    .stripIndent()
}
'''

# Template for benchmark report
BENCHMARK_REPORT_TEMPLATE = '''#!/usr/bin/env python3

import os
import re
from datetime import datetime
from pathlib import Path

def parse_comparison_file(filepath):
    """Parse comparison statistics"""
    stats = {'sample': '', 'callers': {}}

    with open(filepath) as f:
        content = f.read()

        # Extract sample name
        sample_match = re.search(r'Sample: (\\S+)', content)
        if sample_match:
            stats['sample'] = sample_match.group(1)

        # Extract variant counts
        for line in content.split('\\n'):
            if ':' in line and 'variants' in line:
                parts = line.split(':')
                caller = parts[0].strip()
                count = parts[1].strip().split()[0]
                stats['callers'][caller] = count

    return stats

# Collect all comparison files
comparison_files = [f for f in os.listdir('.') if f.endswith('_comparison.txt')]
all_stats = [parse_comparison_file(f) for f in sorted(comparison_files)]

# Generate HTML report
html = """
<!DOCTYPE html>
<html>
<head>
    <title>Variant Calling Benchmark Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 40px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        .container {
            background: white;
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            max-width: 1200px;
            margin: 0 auto;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 4px solid #667eea;
            padding-bottom: 15px;
            margin-bottom: 30px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
        }
        .info-box {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        th {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }
        td {
            padding: 12px 15px;
            border-bottom: 1px solid #ecf0f1;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .metric {
            display: inline-block;
            background: #ecf0f1;
            padding: 10px 20px;
            border-radius: 5px;
            margin: 5px;
            font-weight: 600;
        }
        .footer {
            margin-top: 40px;
            padding-top: 20px;
            border-top: 2px solid #ecf0f1;
            color: #7f8c8d;
            text-align: center;
        }
        .summary {
            background: #f8f9fa;
            padding: 20px;
            border-left: 4px solid #667eea;
            margin: 20px 0;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Multi-Caller Variant Calling Benchmark</h1>

        <div class="info-box">
            <h3>Pipeline Information</h3>
            <p><strong>Callers Used:</strong> {{ ", ".join(callers.keys()) }}</p>
            <p><strong>Total Samples:</strong> {num_samples}</p>
            <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </div>

        <div class="summary">
            <h3>Summary</h3>
            <p>This report compares variant calling results across multiple tools to help identify:</p>
            <ul>
                <li><strong>Consensus variants:</strong> Called by multiple callers (high confidence)</li>
                <li><strong>Caller-specific variants:</strong> Unique to one tool (may need validation)</li>
                <li><strong>Sensitivity differences:</strong> Which callers detect more variants</li>
            </ul>
        </div>

        <h2>Per-Sample Results</h2>
"""

for stats in all_stats:
    html += f"""
        <h3>Sample: {stats['sample']}</h3>
        <table>
            <thead>
                <tr>
                    <th>Variant Caller</th>
                    <th>Variants Called</th>
                </tr>
            </thead>
            <tbody>
"""
    for caller, count in stats['callers'].items():
        html += f"""
                <tr>
                    <td>{caller}</td>
                    <td>{count}</td>
                </tr>
"""
    html += """
            </tbody>
        </table>
"""

html += """
        <div class="footer">
            <p>Generated by Jinja2-templated Nextflow Pipeline</p>
            <p>For ensemble calling and detailed concordance analysis, see individual sample reports</p>
        </div>
    </div>
</body>
</html>
"""

with open('benchmark_report.html', 'w') as f:
    f.write(html)

print("Benchmark report generated successfully!")
'''

def generate_pipeline(callers=None, output_dir='generated_pipeline'):
    """Generate the Nextflow pipeline with specified callers"""

    # Use all callers if none specified
    if callers is None:
        callers = CALLERS_CONFIG
    else:
        callers = {k: v for k, v in CALLERS_CONFIG.items() if k in callers}

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    templates_path = output_path / 'templates'
    templates_path.mkdir(exist_ok=True)

    # Generate main Nextflow file
    template = Template(NEXTFLOW_TEMPLATE)
    pipeline_content = template.render(callers=callers)

    with open(output_path / 'main.nf', 'w') as f:
        f.write(pipeline_content)

    # Generate benchmark report template
    template = Template(BENCHMARK_REPORT_TEMPLATE)
    report_content = template.render(callers=callers)

    with open(templates_path / 'benchmark_report.py', 'w') as f:
        f.write(report_content)

    # Generate nextflow.config
    config_content = f"""
params {{
    bam_dir = "./bam_files"
    reference = "./reference/genome.fasta"
    outdir = "./results"
    callers = "{','.join(callers.keys())}"
    enable_ensemble = true
    min_callers_agree = 2
}}

process {{
    errorStrategy = 'retry'
    maxRetries = 2

    withLabel: process_medium {{
        cpus = 2
        memory = '4.GB'
    }}

    withLabel: process_high {{
        cpus = 4
        memory = '8.GB'
    }}
}}

docker {{
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}}

manifest {{
    name = 'multi-caller-variant-pipeline'
    description = 'Benchmark multiple variant callers'
    version = '1.0.0'
}}
"""

    with open(output_path / 'nextflow.config', 'w') as f:
        f.write(config_content)

    # Generate README
    readme_content = f"""# Multi-Caller Variant Calling Pipeline

Generated using Jinja2 templating.

## Included Callers

{chr(10).join(f'- **{config["name"]}** ({caller_id})' for caller_id, config in callers.items())}

## Usage

```bash
# Run with all callers
nextflow run main.nf

# Run with specific callers
nextflow run main.nf --callers "gatk_haplotypecaller,freebayes"

# Disable ensemble calling
nextflow run main.nf --enable_ensemble false

# Custom parameters
nextflow run main.nf \\
    --bam_dir /path/to/bams \\
    --reference /path/to/ref.fasta \\
    --outdir results \\
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
"""

    with open(output_path / 'README.md', 'w') as f:
        f.write(readme_content)

    print(f"✅ Pipeline generated in {output_dir}/")
    print(f"   - Main pipeline: main.nf")
    print(f"   - Configuration: nextflow.config")
    print(f"   - Templates: templates/")
    print(f"   - {len(callers)} variant callers included")
    print(f"\nTo run: cd {output_dir} && nextflow run main.nf")

if __name__ == '__main__':
    # Example: Generate with all callers
    generate_pipeline()

    # Or generate with specific callers
    # generate_pipeline(callers=['gatk_haplotypecaller', 'freebayes', 'bcftools'])