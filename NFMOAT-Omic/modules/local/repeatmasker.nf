process REPEATMASKER{
    tag "$meta.id"
    label "process_medium"

    container "${ 'https://singularity/tetools' }"

    input:
        tuple val(meta),path(species),path(genome)

    output:
        tuple val(meta),path("*.out"),emit:fasta

    script:
        """
            BuildDatabase -name test  ${genome}

            RepeatModeler -engine ncbi -threads 30 -database test

            cat test-families.fa ${species} > repeat_db.fa
             
            RepeatMasker -xsmall -gff -html -lib repeat_db.fa -pa 30 ${genome} 
        
        """

}
