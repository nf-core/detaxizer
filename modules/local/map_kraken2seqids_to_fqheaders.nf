process MAP_KRAKEN2SEQIDS_TO_FQHEADERS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*kraken2.map.txt'), emit: mapping
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - ${reads} ${meta.id}.kraken2.map.txt <<'PY'
import gzip, sys

def open_file(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')

reads = sys.argv[1:-1]
output = sys.argv[-1]

with open(output, 'w') as out:
    if len(reads) == 2:
        with open_file(reads[0]) as f1, open_file(reads[1]) as f2:
            while True:
                h1 = f1.readline().rstrip()
                h2 = f2.readline().rstrip()
                if not h1 or not h2:
                    break
                f1.readline(); f1.readline(); f1.readline()
                f2.readline(); f2.readline(); f2.readline()
                kid = h1[1:].split()[0].replace('/1','').replace('/2','')
                out.write(f"{kid}\\t{h1[1:]}\\t{h2[1:]}\\n")
    else:
        with open_file(reads[0]) as f1:
            while True:
                h1 = f1.readline().rstrip()
                if not h1:
                    break
                f1.readline(); f1.readline(); f1.readline()
                kid = h1[1:].split()[0].replace('/1','').replace('/2','')
                out.write(f"{kid}\\t{h1[1:]}\\n")
PY

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
