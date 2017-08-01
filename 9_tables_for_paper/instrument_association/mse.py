#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import os
import pandas as pd

def main():

    # Load instruments
    instr = pd.read_csv("../../data/instruments/Myopia/myopia_23andme_variants_170728.tsv",
                        sep="\t", header=0)
    print instr.head()

    # Load associations
    assoc = pd.read_csv("../../5_variant_plots/output/myopiaInstr_stats.tsv",
                        sep="\t", header=0)
    print assoc.head()

    # Merge
    merged = pd.merge(instr, assoc, left_on="MarkerName", right_on="rsid",
                      how="inner")
    print merged.head()

    # Save
    merged.to_csv("myopia_instruments_assocs.tsv", sep="\t", index=None)


    return 0

if __name__ == '__main__':

    main()
