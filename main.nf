params.fasta = "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
params.chrom_sizes = "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes > refs/sacCer3.chrom.sizes"

params.ids = ['SRR3033157', 'SRR3033156'] //, 'SRR3033155'

ch_fastq = channel.fromFilePairs('/scratch/applied-genomics/chipseq/SRR*_{1,2}.fastq')

process create_genome {
    conda 'bwa samtools'
    input:
    path fasta_tar

    output:
    path 'genome.fa', emit: genome_fasta

    """
    tar zxvf ${fasta_tar}
    cat chr*.fa > genome.fa
    """
}

process bwa_index {
    conda 'bwa samtools'
    input:
    path REF 

    output:
    path 'ref', emit: ref
    path '*.fai', emit: fai

    """
    mkdir ref

    bwa index  -p ref/${REF.baseName} $REF
    samtools faidx $REF
    """
}

process bwa_align {
    conda 'bwa samtools'

    input:
    path ref
    tuple val(sample_id), path(reads)

    """"
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    
    bwa mem \$INDEX \\
        $reads \\
        | samtools sort > ${sample_id}.bam
    """
}

process PEAK_ANNOTATION {
    conda 'bwa samtools'

    """
    #!/usr/bin/env/ Rscript

    library(ChIPseeker)
}

workflow {
    create_genome(file(params.fasta))
    bwa_index(create_genome.out.genome_fasta)
    bwa_align(bwa_index.out.ref, ch_fastq)
}
    
