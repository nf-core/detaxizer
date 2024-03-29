process SUMMARY_KRAKEN2 {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10.4 pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
    tuple val(meta),path(kraken2)


    output:
    tuple val(meta), path("*.kraken2_summary.tsv")  , emit: summary
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import subprocess
    import pandas as pd

    def sort_list_of_files_by_pattern(kraken2_complete_list):
        kraken2_dict = {
            "kraken2": [],
            "isolatedkraken2": []
        }
        for entry in kraken2_complete_list:
            if 'classifiedreads.txt' in entry:
                kraken2_dict["kraken2"].append(entry)
            else:
                kraken2_dict["isolatedkraken2"].append(entry)

        return kraken2_dict

    def calculate_lines_of_file(path):
        lines = 0
        with open(path,'r') as f:
            for line in f:
                if line == '\\n':
                    pass
                else:
                    lines += 1
        return lines

    def get_version():
        version_output = subprocess.getoutput('python --version')
        return version_output.split()[1]

    list_files = "${kraken2}".split(" ")
    kraken2_dict = sort_list_of_files_by_pattern(list_files)
    kraken2_dict_lines = {
        "kraken2": 0,
        "isolatedkraken2": 0
    }
    for key in kraken2_dict.keys():
        for entry in kraken2_dict[key]:
            kraken2_dict_lines[key] += calculate_lines_of_file(entry)
        kraken2_dict_lines[key] = [ kraken2_dict_lines[key] ]
    df = pd.DataFrame(kraken2_dict_lines)

    df.index = ["${meta.id}"]

    df.to_csv("${meta.id}.kraken2_summary.tsv",sep='\\t')

    # Generate the version.yaml for MultiQC
    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
