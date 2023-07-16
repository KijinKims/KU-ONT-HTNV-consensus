workflow consensus {
    take:
        bam
        seg
    main:
        consensus_process(bam, seg)
}

process consensus_process {
    conda 'htnv_consensus'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy'
    tag "${params.prefix}:consensus"

    input:
        path bam
        path seg
    output:
        path "${bam.simpleName}.consensus.fasta"
    """
    samtools index $bam
    medaka consensus --model r941_prom_variant_g360 --threads 8 --batch_size 100 $bam ${bam.simpleName}.hdf
    medaka variant $seg ${bam.simpleName}.hdf ${bam.simpleName}.medaka.vcf
    medaka tools annotate --dpsp ${bam.simpleName}.medaka.vcf $seg $bam ${bam.simpleName}.medaka.annotated.vcf
    bcftools view -i 'QUAL>${params.variant_quality_threshold} & INFO/DP>${params.variant_depth_threshold}' ${bam.simpleName}.medaka.annotated.vcf > ${bam.simpleName}.medaka.filtered.vcf
    python ${params.indel_filter_script} --input ${bam.simpleName}.medaka.filtered.vcf --output ${bam.simpleName}.medaka.indel_filtered.vcf
    bgzip -f ${bam.simpleName}.medaka.indel_filtered.vcf
    tabix -p vcf ${bam.simpleName}.medaka.indel_filtered.vcf.gz

    # apply variants to reference to generate consensus
    cat $seg | bcftools consensus ${bam.simpleName}.medaka.indel_filtered.vcf.gz > ${bam.simpleName}.consensus.nodropped.fasta

    # drop low coverage region
    bedtools genomecov -bga -ibam ${bam.simpleName}.bam | awk '\$4 < ${params.low_cov_threshold}' | bedtools merge -i - > low_cov_${params.low_cov_threshold}.bed
    bedtools maskfasta -fi ${bam.simpleName}.consensus.nodropped.fasta -bed low_cov_${params.low_cov_threshold}.bed -fo ${bam.simpleName}.consensus.fasta

    # change fasta header
    header=">${params.prefix}_${seg.simpleName}_consensus low_cov_thrs=${params.low_cov_threshold} var_qual_thrs=${params.variant_quality_threshold} var_depth_thrs=${params.variant_depth_threshold}"
    sed -i "1s/.*/\$header/" ${bam.simpleName}.consensus.fasta
    """
}
