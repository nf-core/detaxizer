process SUMMARIZER {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
        path(tosummarize)
        

    output:
        path("summary.tsv"), emit: summary
        //path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env python
    import glob
    import pandas as pd
    files_kraken2 = glob.glob('*.kraken2_summary.tsv')
    files_blastn = glob.glob('*.blastn_summary.tsv')
    df_kraken2 = pd.DataFrame()
    df_blastn = pd.DataFrame()

    for file in files_kraken2:
        df_local = pd.read_csv(file, sep="\\t", index_col=0)
        df_kraken2 = pd.concat([df_kraken2,df_local])
    for file in files_blastn:
        df_local = pd.read_csv(file, sep="\\t", index_col=0)
        df_blastn = pd.concat([df_blastn,df_local])
    
    df = pd.concat([df_kraken2, df_blastn.reindex(df_kraken2.index)],axis=1)
    df.to_csv("summary.tsv",sep="\\t")
        
    """
}