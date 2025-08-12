process SNPEFF_ANN_VCF{
    tag "$meta.id"
    label 'process_high'

    errorStrategy { task.exitStatus in 150..152 ? 'retry' : 'terminate' }
    maxRetries 10

    container "${ 'https://singularity/snpEFF' }"

    input:
    tuple val(meta),path(smp_vcf),path(fasta),path(gff)

    output:
    tuple val(meta),path("*.ann.vcf.gz"),path("*.txt"),path("*.html"),path("*.results"),emit:ann_results

    script:
    """
    fasta_name=\$(basename $fasta .sort.fa)
    mkdir -p data/abc
    cp $fasta ./data/abc/sequences.fa && cp  $gff ./data/abc/genes.gff 
    java -jar /opt/snpEff/snpEff/snpEff.jar build -c /opt/snpEff/snpEff/snpEff.config  -dataDir \$(pwd)/data -gff3 -v abc -noCheckCds -noCheckProtein 
    java -jar /opt/snpEff/snpEff/snpEff.jar -dataDir \$(pwd)/data -noDownload -c /opt/snpEff/snpEff/snpEff.config abc ${smp_vcf} >\$(basename ${smp_vcf} .vcf).ann.vcf.gz
    perl -alne 'next if \$_ =~ /^#/; \$F[7] =~ /(ANN=\\S+)/; print "\$F[2]\\t\$F[5]\\t\$1"' \$(basename ${smp_vcf} .vcf).ann.vcf.gz >\$(basename ${smp_vcf} .vcf).ann.results
    """
}
