process ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(kraken2results), path(tax2filter)

    output:
    tuple val(meta), path('*classified.txt'), emit: classified
    tuple val(meta), path('*ids.txt')       , emit: classified_ids
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import subprocess
    import re
    from collections import defaultdict

    # taken from compagnon of cogdat https://github.com/ckaipf/compagnon/blob/main/compagnon/service_layer/executions/cogdat/contamination.py
    def parse_kraken_lca_mapping(string: str) -> defaultdict[int, int]:
        d: defaultdict[int, int] = defaultdict(int)
        for kmer in string.strip().split(" "):
            splitted = kmer.strip().split(":")
            if all(x not in ["A", "|"] for x in splitted):
                ncbi_id, kmer_count = (int(x) for x in splitted)
                d[ncbi_id] += kmer_count
        return d

    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    tax2filter = []
    with open('${tax2filter}','r') as file:
        for line in file:
            line = int(line.strip('\\n'))
            tax2filter.append(line)

    filterList = []
    with open('${kraken2results}', 'r') as file:
        with open('${meta.id}.classified.txt', 'w') as outfile:
            for line in file:
                line = line.split("\\t")
                lca_mapping = parse_kraken_lca_mapping(line[4])
                tax2keep_in_this_line = []
                for key in lca_mapping.keys():
                    if key not in tax2filter and key != 0:
                        tax2keep_in_this_line.append(key)
                sum_to_keep = sum([lca_mapping[id_] for id_ in tax2keep_in_this_line])
                sum_to_filter = sum([lca_mapping[id_] for id_ in tax2filter])
                unclassified = lca_mapping[0]

                if (
                    sum_to_filter > $params.cutoff_tax2filter
                    and sum_to_filter/(sum_to_keep + sum_to_filter) > $params.cutoff_tax2keep
                    and sum_to_filter/(unclassified + sum_to_filter) > $params.cutoff_unclassified
                ):
                    filterList.append(line[1])
                    outfile.write("\\t".join(line))

    with open('${meta.id}.ids.txt', 'w') as outfile:
        for entry in filterList:
            outfile.write(entry+"\\n")

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
    """
}
