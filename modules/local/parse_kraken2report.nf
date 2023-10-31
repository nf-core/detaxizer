process PARSE_KRAKEN2REPORT {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"
    
    input:
        tuple val(meta), path(kraken2report)
    
    output:
        path "versions.yml", emit: versions
        path "*.txt", emit: txt
    
    script:
    """
    #!/usr/bin/env python
import subprocess


def read_in_kraken2report(report):
    kraken_list = []

    with open(report, 'r') as f:
        for line in f:
            if not line:
                continue

            line = line.strip('\\n').split('\\t')
            kraken_list.append(line)

    return kraken_list


def kraken_taxonomy2hierarchy(kraken_list):
    tax_dict = {}
    stack = []

    for line in kraken_list:
        level = (len(line[5]) - len(line[5].strip())) / 2
        node = line[5].strip()
        id = line[4]

        while len(stack) > level:
            stack.pop()

        if level == 0 or not stack:
            tax_dict[node] = {"id": id}
            stack = [tax_dict[node]]
        else:
            current_level = stack[-1]
            current_level[node] = {"id": id}
            stack.append(current_level[node])


    return tax_dict


def get_all_keys(nested_dict, parent_key='', sep='\\t'):
    keys = []

    for k, v in nested_dict.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k

        if k == 'id':
            keys.append(parent_key)
        elif isinstance(v, dict):
            keys.extend(get_all_keys(v, new_key, sep=sep))

    return keys


def remove_incomplete_taxa(list_kraken):
    # remove the incomplete paths to the leafs, keep all paths that are ending in a leaf
    list_kraken.sort()
    to_remove = set()

    for i in range(len(list_kraken) - 1):
        if list_kraken[i + 1].startswith(list_kraken[i] + ","):
            to_remove.add(i)

    return [item for idx, item in enumerate(list_kraken) if idx not in to_remove]


def generate_dict_for_lookup(removed_list):
    lookup_dict = {}

    for entry in removed_list:
        new_list = entry.split('\\t')

        while new_list:
            new_entry = new_list.pop(0)

            if new_entry in lookup_dict:
                lookup_dict[new_entry].update(new_list)
            else:
                lookup_dict[new_entry] = set(new_list)

    for key in lookup_dict:
        lookup_dict[key] = sorted(list(lookup_dict[key]))

    sorted_keys = sorted(lookup_dict)
    return {k: lookup_dict[k] for k in sorted_keys}

def get_keys_with_ids(kraken_dict, results=None):
    if results is None:
        results = {}

    for k, v in kraken_dict.items():
        if isinstance(v, dict):
            if 'id' in v:
                results[k] = v['id']
            results.update(get_keys_with_ids(v, results))
    return results

def get_version():
    version_output = subprocess.getoutput('python --version')
    return version_output.split()[1]


kraken_list = read_in_kraken2report('$kraken2report')
kraken_dict = kraken_taxonomy2hierarchy(kraken_list)
list_kraken = get_all_keys(kraken_dict)
result = remove_incomplete_taxa(list_kraken)
result = generate_dict_for_lookup(result)
result_to_filter = result["$params.tax2filter"] + ["$params.tax2filter"]
result = get_keys_with_ids(kraken_dict)
result_tax_name_id = []
for entry in result_to_filter:
    tax_name = entry
    tax_id = str(result[tax_name])
    #tax_name_id = tax_name + " (taxid " + tax_id + ")"
    result_tax_name_id.append(tax_id)
    
with open('taxa_to_filter.txt', "w") as f:
    for entry in result_tax_name_id:
        f.write(entry+'\\n')

# Generate the version.yaml for MultiQC
with open('versions.yml', 'w') as f:
    f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
    f.write(f'    python: {get_version()}\\n')

    """
}