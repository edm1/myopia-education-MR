#!/bin/sh
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:01:00:00
#PBS -N n-AS-ukbbvis


# Change to work dir if set
if [ ! -z ${PBS_O_WORKDIR+x} ]; then
	cd $PBS_O_WORKDIR
fi

# module load apps/qctool-2.0
module load apps/qctool-1.4

#
# Set variables
#

# Working directory
wd=/panfs/panasas01/sscm/em13383/phd/ukbiobank/myopia_EA_project
# Genotype pattern
arg_geno_bgen=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_imp_chrCHROM_v2.bgen
arg_geno_bgen_index=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/ukb_bgi_chrCHROM_v2.bgi
# Sample file
arg_sample_file=$wd/data/ukb878_imp_chr1_v2_s487406.sample
# EA instruments
arg_ea_instruments=$wd/data/instruments/Educational_attainment/EA_Okbay_variants_inHRC_170728.tsv
# Myopia instruments
arg_myopia_instruments=$wd/data/instruments/Myopia/myopia_23andme_variants_inHRC_170728.tsv
# Refractive error instruments
arg_re_instruments=$wd/data/instruments/Refractive_error/cream_variants_2.tsv
# Visual acuity instruments
arg_va_instruments=$wd/data/instruments/Visual_acuity/VA_ALSPAC.tsv


#
# Generate allele scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

mkdir -p output

#
# Do EA
#

outpref=temp_EA

# Get SNP list
snp_list="$outpref"_snplist.txt
awk 'NR != 1 {print $1}' $arg_ea_instruments | sort > $snp_list

# Extract
for chrom in {1..22}; do
    inbgen=${arg_geno_bgen/CHROM/$chrom}
    inbgenidx=${arg_geno_bgen_index/CHROM/$chrom}
    bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list -v11 > $outpref.$chrom.bgen
done

# Convert to gen
inbgen=$outpref.#.bgen
qctool -g $inbgen -maf 0.01 1 -info 0.3 1 -og - | gzip -c > $outpref.gen.gz

# Identify missing variants
zcat $outpref.gen.gz | cut -f 3 -d " " | sort | comm -23 $snp_list - > $outpref.missingVariants.txt

# Convert gen format to alleleB dosage
zcat $outpref.gen.gz | pypy gen_to_expected.py > $outpref.dosages

# Merge weights with dosages to produce final output
outname=output/EA_UKBB_Okbay.dosage
python merge_weights_to_dosage.py $outpref.dosages $arg_ea_instruments \
    $arg_sample_file $outname
