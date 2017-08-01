#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import pandas as pd
from functools import partial
from pprint import pprint

def main():

    # Args
    indosage  = sys.argv[1]
    inweights = sys.argv[2]
    insample = sys.argv[3]
    outf = sys.argv[4]

    # Load weights into dict
    weight_dict = load_weights(inweights)
    # pprint(weight_dict)
    # sys.exit()

    # Load sample names
    samples = load_sample_list(insample)

    # Load dosage data
    dosage = pd.read_csv(indosage, sep=" ", header=None)
    dosage.columns = ["chr", "snpid", "rsid", "pos", "alleleA", "alleleB"] + samples

    # Flip weights so that they correspond to 1 increase in alleleB
    weights = dosage.apply(partial(match_allele_weights, w=weight_dict), axis=1)
    dosage.insert(6, "weightB", weights)

    # Write file
    dosage.to_csv(outf, sep="\t", index=None)

    return 0

def load_sample_list(inf):
    """ Loads list of sample names from SNPTEST sample file
    """
    sample_list = []
    with open(inf, "r") as in_h:
        # Skip 2 headers
        in_h.readline()
        in_h.readline()
        # Read all in
        for line in in_h:
            id1, id2, missing = line.rstrip().split(" ")
            if not id1 == id2:
                sys.exit("ID1 and ID2 were different in sample file:", id1, id2)
            sample_list.append(id1)
    return sample_list

def match_allele_weights(row, w):
    """ Check that effect allele matches alleleB, if not flip direction of the
        weight.
    """
    # Get info from dosage row
    rsid = row["rsid"]
    alleleA = row["alleleA"]
    alleleB = row["alleleB"]
    # Check that alleles are the same
    if not sorted([alleleA, alleleB]) == sorted(list(w[rsid]["alleles"])):
        sys.exit(("Error: Alleles don't match for: ", rsid))
    # If effect allele is alleleB, return weight otherwise *-1
    if alleleB == w[rsid]["alleles"][1]:
        return float(w[rsid]["weight"])
    else:
        return -1*float(w[rsid]["weight"])

def load_weights(inf, sep="\t"):
    """ For each variant load weight and (Other_allele, Effect_allele) into
        a dict. Cols required: MarkerName, Beta, Other_allele, Effect_allele
    """
    w = {}
    with open(inf, "r") as in_h:
        # Load header and discover column indexs
        header = in_h.readline().rstrip().split(sep)
        cols = {}
        for key in ["MarkerName", "Beta", "Other_allele", "Effect_allele"]:
            cols[key] = header.index(key)
        # For each line save variants, alleles and weight
        for line in in_h:
            parts = line.rstrip().split(sep)
            alleles = (parts[cols["Other_allele"]], parts[cols["Effect_allele"]])
            w[parts[cols["MarkerName"]]] = {"weight":parts[cols["Beta"]],
                                            "alleles":alleles}
    return w

if __name__ == '__main__':

    main()
