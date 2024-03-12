#! /usr/bin/env python

# Written by Jannik Seidel and released under MIT license.

import argparse
import sys

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Parsing the kraken2 report to get the taxon/taxonomic subtree which should be assessed/filtered.")

    parser.add_argument("-i", "--input", help="Path to input file.", type=str, required=True)

    parser.add_argument("-t", "--tax2filter", help="Pass the taxon to assess/to filter to this flag.", type=str, required=True)

    return parser.parse_args()

def read_in_kraken2report(report):
    krakenList = []
    with open(report, 'r') as f:
        for line in f:
            if not line:
                continue
            lineToAppend = line.replace('\n','').split('\t')
            krakenList.append(lineToAppend)
    return krakenList

def krakenTaxonomy2hierarchy(krakenList):
    # function to read in the taxonomic information from a kraken2report into a nested dictionary
    # init of nested dictonary to store taxonomy
    taxDict = {}
    # init of stack to keep track of the different taxonomic levels of each entry
    stack = []

    # processing of the lines in the krakenList from the kraken2 report
    for line in krakenList:
        # level is indicated by the leading whitespaces
        # (2 per hierarchy level) of the fifth entry in
        # the kraken2 report. To address this to get the
        # right level of the current node the ammount of
        # whitespaces is divided by 2
        level = (len(line[5]) - len(line[5].strip())) / 2
        node = line[5].strip()

        # Taxonomic ID, e.g. '9606' for 'Homo sapiens'
        id = line[4]

        # while the length of the stack is greater than the level of the current node,
        # remove the newest entry of the stack to go to the level of the parent
        # node of the current node. When this is reached, continue
        while len(stack) > level:
            stack.pop()

        # for the initial levels (e.g. '0' for 'unclassified' or '1' for 'root')
        # or if the stack is empty add a new parent node to the phylogenetic tree
        if level == 0 or not stack:
            taxDict[node] = {"id": id}
            stack = [taxDict[node]]
        # Add the subtree to the current level of the tree (inside of the taxDict
        # nested dictionary) and append this part of the taxDict to the stack
        else:
            currentLevel = stack[-1]
            currentLevel[node] = {"id": id}
            stack.append(currentLevel[node])

    return taxDict

def getAllKeys(nestedDict, parentKey='', sep='\t'):
# get all the keys in the nested taxonomic dictionary recursively.
# The format of the returned list is as follows:
# [
#       'unclassified',
#       'root',
#       'root\tcellular organism',
#       'root\tcellular organism\tBacteria',
#       ...
#       'root\tcellular organism\tEukaryota',
#       ...
# ]
    keys = []

    for k, v in nestedDict.items():
        newKey = f"{parentKey}{sep}{k}" if parentKey else k

        if k == 'id':
            keys.append(parentKey)
        elif isinstance(v, dict):
            keys.extend(getAllKeys(v, newKey, sep=sep))

    return keys

def removeIncompleteTaxa(listKrakenTaxonomicKeys, sep = '\t'):
    # remove the incomplete paths to the leafs of the hierarchical tree, keep all paths that are ending in a leaf
    listKrakenTaxonomicKeys.sort()
    to_remove = set()

    for i in range(len(listKrakenTaxonomicKeys) - 1):
        if listKrakenTaxonomicKeys[i + 1].startswith(listKrakenTaxonomicKeys[i] + sep):
            to_remove.add(i)

    return [item for idx, item in enumerate(listKrakenTaxonomicKeys) if idx not in to_remove]

def generateDictForLookupOfTaxonomicSubentries(removedIncompleteTaxaList):
    # This function generates a lookup dictionary that allows to get all subentries
    # of a certain taxonomic entry (e.g. for 'Homo' it looks as follows:
    # {'Homo':'Homo sapiens'})
    lookupDict = {}

    for entry in removedIncompleteTaxaList:
        listOfTaxonomicNames = entry.split('\t')

        while listOfTaxonomicNames:
            newTaxomomicName = listOfTaxonomicNames.pop(0)

            if newTaxomomicName in lookupDict:
                lookupDict[newTaxomomicName].update(listOfTaxonomicNames)
            else:
                lookupDict[newTaxomomicName] = set(listOfTaxonomicNames)

    for key in lookupDict:
        lookupDict[key] = sorted(list(lookupDict[key]))

    sortedKeys = sorted(lookupDict)
    return {k: lookupDict[k] for k in sortedKeys}

def getKeysWithIDs(krakenHierarchicalTaxonomyDict, results=None):
    # itterate over the hierarchical taxonomy dictionary and create
    # a mapping from taxonomic name to taxonomic id.
    # This is done recursively and results=None is the inital state.
    if results is None:
        results = {}

    for k, v in krakenHierarchicalTaxonomyDict.items():
        if isinstance(v, dict):
            if 'id' in v:
                results[k] = v['id']
            results.update(getKeysWithIDs(v, results))
    return results

def main():

    args = parse_args()

    krakenList = read_in_kraken2report(args.input)
    krakenHierarchicalTaxonomyDict = krakenTaxonomy2hierarchy(krakenList)
    listKrakenTaxonomicKeys = getAllKeys(krakenHierarchicalTaxonomyDict)
    resultRemovedIncompleteTaxa = removeIncompleteTaxa(listKrakenTaxonomicKeys)
    resultLookUpDict = generateDictForLookupOfTaxonomicSubentries(resultRemovedIncompleteTaxa)
    try:
        # Try if the taxon to assess for/to filter is in the database (represented by the lookup dictionary)
        # and add all subentries of the 'tax2filter' plus the 'tax2filter' to a single list
        resultToFilter = resultLookUpDict[args.tax2filter] + [args.tax2filter]
    except KeyError:
        raise KeyError("The taxaomic group/taxon you want to check for/filter out is not in the kraken database. Use a database that includes the taxonomic group or taxon or change the tax2filter parameter to something that is in the database.")

    resultLookUpOfIDsForTaxonomicNameDict = getKeysWithIDs(krakenHierarchicalTaxonomyDict)
    resultTaxNameID = []
    for entry in resultToFilter:
        tax_name = entry
        tax_id = str(resultLookUpOfIDsForTaxonomicNameDict[tax_name])
        resultTaxNameID.append(tax_id)

    with open('taxa_to_filter.txt', "w") as f:
        for entry in resultTaxNameID:
            f.write(entry+'\n')

if __name__ == "__main__":
    sys.exit(main())
