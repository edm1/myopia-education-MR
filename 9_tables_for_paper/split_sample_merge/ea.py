#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import os
import pandas as pd

def main():

    # Load sample A
    a = pd.read_csv("../../8_MRBase/using_split_sample/1_split_sample_weights/output/eaInstr_stats_splitA.tsv",
                        sep="\t", header=0)
    print a.head()

    # Load sample B
    b = pd.read_csv("../../8_MRBase/using_split_sample/1_split_sample_weights/output/eaInstr_stats_splitB.tsv",
                        sep="\t", header=0)
    print b.head()

    # Merge
    merged = pd.merge(a, b, on="rsid", how="inner")
    print merged.head()

    # Save
    merged.to_csv("ea_splitSample_merge.tsv", sep="\t", index=None)


    return 0

if __name__ == '__main__':

    main()
