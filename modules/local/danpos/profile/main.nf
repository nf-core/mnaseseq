
process DANPOS_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'docker.io/biosseg/danpos:v1' }"

    input:
        tuple val(meta), path(wig)
        tuple val(meta), path(bed)
        tuple val(meta), path(genefile)

    output:
        tuple val(meta), path("*.bam")  , emit: xls
        tuple val(meta), path("*.pdf")
        path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Creating the gene file from the gtf 
    danpos.py profile \\
        ${wig} \\
        $args \\
        --bed3file_paths ${bed} \\
        --genefile_paths ${genefile} \\
    cat <<-VERSIONS >> versions.yml
    ${task.process}: danpos_profile, danpos \$(danpos.py --version | head -n1)
    VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.pdf
    touch ${prefix}.pdf
    touch version.yml
    """
}
