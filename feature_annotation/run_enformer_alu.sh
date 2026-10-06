#!/bin/bash

# Paths - edit for your environment
SUPREMO_DIR="${SUPREMO_DIR:-/pollard/home/szhang20/akita_variant_scoring}"
AKITA_DIR="${AKITA_DIR:-/pollard/home/szhang20/akita}"
PROJECT_DIR="${ALU_PROJECT_DIR:-/pollard/home/szhang20/alu}"



#to run enformer
#using supremo_enformer_env


#set virtual memory limits to 300gb
ulimit -v $((100 * 1024 * 1024)) 

input_file="${PROJECT_DIR}/data/20231122_hg38Alus_100sample_repName_DEL.txt"
#input_file="${PROJECT_DIR}/data/20250505_10alutest_DEL.txt"

out_dir="${PROJECT_DIR}/results/20250929_alu5100_enformer_tracks/"
out_name="20231122_hg38Alus_100sample_repName_DEL_ENFORMER"
# tracks=" "DNASE:epidermal melanocyte" "DNASE:HFF-Myc originated from foreskin fibroblast" "DNASE:lower leg skin female adult (53 years)" " 


python /pollard/home/szhang20/SuPreMo/scripts/streamlined_SuPreMo_multitracks.py $input_file --dir $out_dir --file $out_name --fa ${SUPREMO_DIR}/data/hg38.fa --nrows 500 --get_Enformer_scores --get_tracks --selected_tracks "CHIP:CTCF:GM12878" "CHIP:CTCF:H1-hESC" "CHIP:CTCF:IMR-90" "CHIP:CTCF:HCT116" "CHIP:CTCF:foreskin fibroblast male newborn" "CHIP:CTCF:astrocyte" "CHIP:CTCF:myotube" "CHIP:CTCF:heart left ventricle female adult (53 years)" "CHIP:CTCF:neural progenitor cell originated from H9" "CHIP:CTCF:astrocyte of the cerebellum" "CHIP:CTCF:gastrocnemius medialis female adult (53 years)" "CHIP:CTCF:transverse colon male adult (37 years)" "CHIP:CTCF:breast epithelium female adult (53 years)" "CHIP:CTCF:adrenal gland male adult (54 years)" "CHIP:CTCF:ovary female adult (51 year)" "CHIP:H3K27me3:GM12878" "CHIP:H3K27me3:H1-hESC" "CHIP:H3K27me3:IMR-90" "CHIP:H3K27me3:HCT116" "CHIP:H3K27me3:foreskin fibroblast male newborn" "CHIP:H3K27me3:neural progenitor cell" "CHIP:H3K27me3:gastrocnemius medialis female adult (53 years)" "CHIP:H3K27me3:transverse colon male adult (37 years)" "CHIP:H3K27me3:breast epithelium female adult (53 years)" "CHIP:H3K27me3:adrenal gland male adult (34 years)" "CHIP:H3K27me3:ovary female adult (30 years)" "DNASE:GM12878" "DNASE:H1-hESC" "DNASE:IMR-90" "DNASE:HCT116" "DNASE:frontal cortex female adult (67 years)" "DNASE:cerebellum male adult (27 years) and male adult (35 years)" "DNASE:frontal cortex male adult (27 years) and male adult (35 years)" "DNASE:foreskin fibroblast male newborn" "DNASE:astrocyte" "DNASE:myotube originated from skeletal" "DNASE:neural progenitor cell originated from H9" "DNASE:heart left ventricle female adult (53 years)" "DNASE:gastrocnemius medialis female adult (51 year)" "DNASE:transverse colon male adult (54 years)" "DNASE:breast epithelium female adult (51 year)" "DNASE:adrenal gland male adult (54 years)" "DNASE:ovary female adult (51 year)" "CHIP:REST:H1-hESC" "CHIP:REST:GM12878" "CHIP:YY1:GM12878" "CHIP:RAD21:H1-hESC" "CHIP:RAD21:GM12878" "CHIP:RAD21:IMR-90" "CHIP:ZNF143:GM12878" "CHIP:ZNF143:H1-hESC" "CAGE:brain, adult, pool1" "CAGE:colon, adult, pool1" "CAGE:heart, adult, pool1" "CAGE:skeletal muscle, adult, pool1" "CAGE:temporal lobe, adult, pool1" "CAGE:frontal lobe, adult, pool1" "CAGE:heart, fetal, pool1" "CAGE:B lymphoblastoid cell line: GM12878"
