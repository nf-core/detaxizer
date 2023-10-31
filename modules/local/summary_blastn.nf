process SUMMARY_BLASTN {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
        tuple val(meta), path(blastn_1), path(blastn_2), path(blastn_3), path(filteredblastn_1), path(filteredblastn_2), path(filteredblastn_3)
        

    output:
        tuple val(meta), path("*.blastn_summary.tsv"), emit: summary
        //path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    def get_lines_in_file(filename):
        line_counter = 0
        with open(filename,'r') as f:
            for line in f:
                if line == '\\n':
                    continue
                else:
                    line_counter += 1
        return line_counter

    def get_unique_read_ids(filename):
        read_id_set = set()
        with open(filename, 'r') as f:
            for line in f:
                if line == '\\n':
                    continue
                else:
                    line = line.split("\\t")[0].strip("/1").strip("/2")
                    read_id_set.add(line)
        return read_id_set

    def sort_blastn_and_filteredblastn(blastn,filteredblastn):
        if blastn.count("NO") != filteredblastn.count("NO"):
            raise ValueError("Ammount of 'NO' entries do not match in blastn and filteredblastn lists.")
        else:
            blastnsummary_dict = {}
            if "NO" in blastn and blastn.count("NO") != 3:
                blastn.remove("NO")
                blastn.sort()
                filteredblastn.remove("NO")
                filteredblastn.sort()
            
            for entry in blastn:
                if "_R1" in entry:
                    blastnsummary_dict["blastn_1"] = entry
                elif "_R2" in entry:
                    blastnsummary_dict["blastn_2"] = entry
                elif "_R3" in entry:
                    blastnsummary_dict["blastn_3"] = entry
            
            for entry in filteredblastn:
                if "_R1" in entry:
                    blastnsummary_dict["filteredblastn_1"] = entry
                elif "_R2" in entry:
                    blastnsummary_dict["filteredblastn_2"] = entry
                elif "_R3" in entry:
                    blastnsummary_dict["filteredblastn_3"] = entry
            
            summary_keys = blastnsummary_dict.keys()
            if "blastn_1" not in summary_keys and "filteredblastn_1" not in summary_keys:
                blastnsummary_dict["blastn_1"] = "NA"
                blastnsummary_dict["filteredblastn_1"] = "NA"
            if "blastn_2" not in summary_keys and "filteredblastn_2" not in summary_keys:
                blastnsummary_dict["blastn_2"] = "NA"
                blastnsummary_dict["filteredblastn_2"] = "NA"
            if "blastn_3" not in summary_keys and "filteredblastn_3" not in summary_keys:
                blastnsummary_dict["blastn_3"] = "NA"
                blastnsummary_dict["filteredblastn_3"] = "NA"
            return blastnsummary_dict



    blastn = [
        "${blastn_1}".strip("_FILE1").strip("_FILE2").strip("_FILE3"),
        "${blastn_2}".strip("_FILE1").strip("_FILE2").strip("_FILE3"),
        "${blastn_3}".strip("_FILE1").strip("_FILE2").strip("_FILE3")
    ]

    filteredblastn = [
        "${filteredblastn_1}".strip("_FILE4").strip("_FILE5").strip("_FILE6"),
        "${filteredblastn_2}".strip("_FILE4").strip("_FILE5").strip("_FILE6"),
        "${filteredblastn_3}".strip("_FILE4").strip("_FILE5").strip("_FILE6")
    ]

    blastnsummary_dict = sort_blastn_and_filteredblastn(blastn,filteredblastn)
    unique_ids_blastn = set()
    unique_ids_filteredblastn = set()
    lines_blastn = 0
    lines_filtered_blastn = 0
    counter_NA = 0
    for key in blastnsummary_dict.keys():
        if blastnsummary_dict[key] == "NA":
            blastnsummary_dict[key] = pd.NA
            counter_NA += 1
        elif "filteredblastn" in key:
            lines_filtered_blastn += get_lines_in_file(blastnsummary_dict[key])
            unique_ids_filteredblastn.update(get_unique_read_ids(blastnsummary_dict[key]))
        else:
            lines_blastn += get_lines_in_file(blastnsummary_dict[key])
            unique_ids_blastn.update(get_unique_read_ids(blastnsummary_dict[key]))
    
    final_summary_blastnfilteredblastn = {}
    if counter_NA == 6:
        final_summary_blastnfilteredblastn["blastn_unique_ids"] = pd.NA
        final_summary_blastnfilteredblastn["blastn_lines"] = pd.NA
        final_summary_blastnfilteredblastn["filteredblastn_unique_ids"] = pd.NA
        final_summary_blastnfilteredblastn["filteredblastn_lines"] = pd.NA
    else:
        final_summary_blastnfilteredblastn["blastn_unique_ids"] = len(unique_ids_blastn)
        final_summary_blastnfilteredblastn["blastn_lines"] = lines_blastn
        final_summary_blastnfilteredblastn["filteredblastn_unique_ids"] = len(unique_ids_filteredblastn)
        final_summary_blastnfilteredblastn["filteredblastn_lines"] = lines_filtered_blastn

    df = pd.DataFrame(final_summary_blastnfilteredblastn, index = ["${meta.id}"])

    df.to_csv("${meta.id}.blastn_summary.tsv", sep="\\t")

    """
    
    
}