process SUMMARIZER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10.4 pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
    tuple val(meta), path(tosummarize)

    output:
    tuple val(meta), path("summary.tsv"), emit: summary
    path("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python
    import glob
    import pandas as pd
    import subprocess
    import numpy

    def get_version():
        version_output = subprocess.getoutput('python --version')
        return version_output.split()[1]

    files_kraken2 = glob.glob('*.kraken2_summary.tsv')
    files_blastn = glob.glob('*.blastn_summary.tsv')

    kraken2_dfs = [pd.read_csv(file, sep="\\t", index_col=0) for file in files_kraken2]
    df_kraken2 = pd.concat(kraken2_dfs)
    if files_blastn != []:
        blastn_dfs = [pd.read_csv(file, sep="\\t", index_col=0) for file in files_blastn]
        df_blastn = pd.concat(blastn_dfs)
        summary_df = df_kraken2.join(df_blastn)
        summary_df.to_csv("summary.tsv", sep="\\t")
    else:
        summary_df = df_kraken2
        summary_df.to_csv("summary.tsv", sep="\\t")

    # Generate the version.yaml for MultiQC
    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
        f.write(f'    numpy: {numpy.__version__}\\n')
    """
}
