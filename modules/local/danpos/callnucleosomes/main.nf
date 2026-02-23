

process DANPOS_CALLNUCLEOSOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' : 
    'docker.io/biosseg/danpos:v1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.xls"), emit: xls
    tuple val(meta), path("*.wig"), emit: wig
    path 'versions.yml',             emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    danpos.py dpos \\
        ${bam} \\
        --out ./result/ \\
        $args > out.txt

    mv ./result/*/* .

    cat <<-VERSIONS >> versions.yml
    ${task.process}: danpos_callnucleosome, danpos \$(danpos.py --version | head -n1)
    VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.xls
    touch ${prefix}.wig
    touch versions.yml
    """
}
