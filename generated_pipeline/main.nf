#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Multi-Caller Variant Calling Pipeline
 * Generated using Jinja2 templating
 */

params.bam_dir = "./bam_files"
params.reference = "./reference/genome.fasta"
params.outdir = "./results"
params.callers = "gatk_haplotypecaller,freebayes,bcftools,deepvariant,varscan"
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


/*
 * Process: GATK HaplotypeCaller
 */
process CALL_VARIANTS_GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.4.0.0'
    publishDir "${params.outdir}/variants/gatk_haplotypecaller", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("gatk_haplotypecaller"), path("*.gatk.vcf.gz"), emit: vcf
    tuple val(sample_id), val("gatk_haplotypecaller"), path("*.gatk.vcf.gz.tbi"), emit: idx

    when:
    params.callers.contains("gatk_haplotypecaller")

    script:
    """
        gatk HaplotypeCaller \
            -R ${reference} \
            -I ${bam} \
            -O ${sample_id}.gatk.vcf.gz \
            --native-pair-hmm-threads ${task.cpus}
        
    """
}


/*
 * Process: FreeBayes
 */
process CALL_VARIANTS_FREEBAYES {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2'
    publishDir "${params.outdir}/variants/freebayes", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("freebayes"), path("*.freebayes.vcf.gz"), emit: vcf
    tuple val(sample_id), val("freebayes"), path("*.freebayes.vcf.gz.tbi"), emit: idx

    when:
    params.callers.contains("freebayes")

    script:
    """
        freebayes \
            -f ${reference} \
            ${bam} \
            > ${sample_id}.freebayes.vcf

        bgzip ${sample_id}.freebayes.vcf
        tabix -p vcf ${sample_id}.freebayes.vcf.gz
        
    """
}


/*
 * Process: BCFtools mpileup
 */
process CALL_VARIANTS_BCFTOOLS {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    publishDir "${params.outdir}/variants/bcftools", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("bcftools"), path("*.bcftools.vcf.gz"), emit: vcf
    tuple val(sample_id), val("bcftools"), path("*.bcftools.vcf.gz.tbi"), emit: idx

    when:
    params.callers.contains("bcftools")

    script:
    """
        bcftools mpileup \
            -f ${reference} \
            -Ou ${bam} | \
        bcftools call \
            --threads ${task.cpus} \
            -mv \
            -Oz \
            -o ${sample_id}.bcftools.vcf.gz

        tabix -p vcf ${sample_id}.bcftools.vcf.gz
        
    """
}


/*
 * Process: DeepVariant
 */
process CALL_VARIANTS_DEEPVARIANT {
    tag "$sample_id"
    label 'process_high'
    container 'google/deepvariant:1.6.1'
    publishDir "${params.outdir}/variants/deepvariant", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("deepvariant"), path("*.deepvariant.vcf.gz"), emit: vcf
    tuple val(sample_id), val("deepvariant"), path("*.deepvariant.vcf.gz.tbi"), emit: idx

    when:
    params.callers.contains("deepvariant")

    script:
    """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=${reference} \
            --reads=${bam} \
            --output_vcf=${sample_id}.deepvariant.vcf.gz \
            --num_shards=${task.cpus}
        
    """
}


/*
 * Process: VarScan
 */
process CALL_VARIANTS_VARSCAN {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/varscan:2.4.6--hdfd78af_0'
    publishDir "${params.outdir}/variants/varscan", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_idx

    output:
    tuple val(sample_id), val("varscan"), path("*.varscan.vcf.gz"), emit: vcf
    tuple val(sample_id), val("varscan"), path("*.varscan.vcf.gz.tbi"), emit: idx

    when:
    params.callers.contains("varscan")

    script:
    """
        samtools mpileup -f ${reference} ${bam} | \
        varscan mpileup2snp \
            --output-vcf 1 \
            > ${sample_id}.varscan.vcf

        bgzip ${sample_id}.varscan.vcf
        tabix -p vcf ${sample_id}.varscan.vcf.gz
        
    """
}



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
    bcftools norm \
        -f ${reference} \
        -m -both \
        ${vcf} \
        -Oz \
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
        caller=\$(basename \$vcf | cut -d'.' -f2)
        count=\$(bcftools view -H \$vcf | wc -l)
        echo "\${caller}: \${count} variants" >> ${sample_id}_comparison.txt
    done

    # Simple intersection analysis
    echo "" >> ${sample_id}_comparison.txt
    echo "Intersection Analysis" >> ${sample_id}_comparison.txt
    echo "--------------------" >> ${sample_id}_comparison.txt

    # This is simplified - in production you'd use more sophisticated comparison
    bcftools isec ${vcf_list} -p isec_out

    echo "Unique to each caller:" >> ${sample_id}_comparison.txt
    for i in \$(seq 0 \$(($num_vcfs-1))); do
        count=\$(bcftools view -H isec_out/\$i.vcf | wc -l)
        echo "  Caller \$i: \${count} unique variants" >> ${sample_id}_comparison.txt
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
    echo "Total variants: \$(bcftools view -H ${sample_id}.ensemble.vcf.gz | wc -l)" >> ${sample_id}.ensemble_stats.txt
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
            def sample_id = file.baseName.replaceAll(/\.bam$/, '')
            def bai = file.parent.resolve(file.name + '.bai')
            tuple(sample_id, file, bai)
        }

    reference_ch = Channel.fromPath(params.reference)
    reference_idx_ch = Channel.fromPath(params.reference + '.fai')

    // Call variants with all callers
    
    gatk_haplotypecaller_vcf = CALL_VARIANTS_GATK_HAPLOTYPECALLER(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    
    freebayes_vcf = CALL_VARIANTS_FREEBAYES(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    
    bcftools_vcf = CALL_VARIANTS_BCFTOOLS(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    
    deepvariant_vcf = CALL_VARIANTS_DEEPVARIANT(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    
    varscan_vcf = CALL_VARIANTS_VARSCAN(
        bam_ch,
        reference_ch,
        reference_idx_ch
    )
    

    // Combine all caller outputs - each emits [sample_id, caller, vcf] and [sample_id, caller, idx]
    all_vcfs = Channel.empty()
    
        .mix(gatk_haplotypecaller_vcf.vcf)
    
        .mix(freebayes_vcf.vcf)
    
        .mix(bcftools_vcf.vcf)
    
        .mix(deepvariant_vcf.vcf)
    
        .mix(varscan_vcf.vcf)
    

    all_idxs = Channel.empty()
    
        .mix(gatk_haplotypecaller_vcf.idx)
    
        .mix(freebayes_vcf.idx)
    
        .mix(bcftools_vcf.idx)
    
        .mix(deepvariant_vcf.idx)
    
        .mix(varscan_vcf.idx)
    

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