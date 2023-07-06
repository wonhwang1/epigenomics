params.fasta = "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
params.chrom_sizes = "http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes > refs/sacCer3.chrom.sizes"

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
    """
    bwa index $REF
    samtools faidx $REF
    """
}

workflow {
    create_genome(file(params.fasta))
    bwa_index(create_genome.out.genome_fasta)
}
    
