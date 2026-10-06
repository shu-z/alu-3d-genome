"""Converted from Alus_feature_overlap.ipynb; logic unchanged."""

import gzip
import math
import os
import re
import subprocess

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import pybedtools

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

#score path
res_dir=f'{PROJECT_DIR}/results/'
score_dir=f'{res_dir}20240308/'
data_dir=f'{PROJECT_DIR}/data/20240308_Alus_sampling/'

# #first look and combine the sampled Alu df and scores 
Alu_sample_path=f'{data_dir}20240308_hg38Alus_1000_random_alignment_sample_DEL.txt'
Alu_sampled=pd.read_csv(Alu_sample_path, sep='\t')
Alu_sampled

rowsize = 2500
score_df_list=[]

nrow_total=len(Alu_sampled)

for i in range(0, nrow_total, rowsize):
#for i in [0]:
    #print('i is ', i)
    df_score=pd.read_csv(f'{score_path}Alu_sample_{i}_scores', sep='\t')
    df_alu=pd.read_csv(f'{data_path}split_files/Alu_sample_{i}.txt', sep='\t')

    #df_score = df_score.drop('var_index', axis=1)

    print(df_score)
    print(df_alu)
    #assert len(df_score) == len(df_alu)
    
    #df_i=pd.concat([df_alu, df_score], axis=1)
    #score_df_list.append(df_i)
    
    #print(len(df_i))

score_df=pd.concat(score_df_list)
score_df=score_df.reset_index(drop=True)
score_df

random_sample_idx=pd.read_csv(f'{data_dir}random_sample_indices.txt', header=None)
random_sample_idx
random_alu_df=score_df[score_df['row_idx'].isin(random_sample_idx[0])]
#read in random alu file 
#random_alu_df=pd.read_csv(f'{score_dir}/combined/20240308_Alu_1000randomsample_scores.txt', sep='\t', index_col=0)
random_alu_df

random_alu_df[random_alu_df['POS']==136537598]

random_alu_df_sample=random_alu_df.sample(n=100, replace=False, random_state=1)
random_alu_df_sample

random_alu_df_sample.to_csv(f'{PROJECT_DIR}/results/20250604_Alu_1000randomsample_100SAMPLECHECK.txt', sep='\t', index=False)

#pick this random state so there's a mix of strands 
random_alu_df_sample5=random_alu_df.sample(n=5, replace=False, random_state=3)
random_alu_df_sample5

random_alu_df_sample5.to_csv(f'{PROJECT_DIR}/results/20250604_Alu_1000randomsample_5SAMPLECHECK.txt', sep='\t', index=False)

random_alu_df_sample5

random_alu_df_sample=pd.read_csv(f'{res_dir}20250604_Alu_1000randomsample_100SAMPLECHECK.txt', sep='\t')

random_sample_scores=pd.read_csv(f'{res_dir}20250604_Alu_1000randomsample_100SAMPLECHECK_scores', sep='\t')
random_sample_scores.columns=random_sample_scores.columns + '_new'
random_sample_all=pd.concat([random_alu_df_sample, random_sample_scores], axis=1)
random_sample_all

alus_5100=pd.read_csv(f'{res_dir}202520513_hg38Alus_100sample_repName_DEL_combinedscores.txt', sep='\t', index_col=0)
alus_5100

merged_df = pd.merge(alus_5100, random_sample_all, on=['CHROM', 'POS', 'END'], how='inner') 
merged_df

plt.scatter(merged_df['mse_mean_x'], merged_df['mse_mean_new'])

supremo_latest_scores=pd.read_csv(f'{res_dir}20250604_Alu_1000randomsample_100SAMPLECHECK_SUPREMO_LATEST_scores', sep='\t')
alu_supremo_scores=pd.read_csv(f'{res_dir}20250604_Alu_1000randomsample_100SAMPLECHECK_alu_supremo_scores', sep='\t')
alu_supremo_mutate_scores=pd.read_csv(f'{res_dir}20250604_Alu_1000randomsample_100SAMPLECHECK_alu_supremo_MUTATE_scores', sep='\t')

#OK THE ISSUE IS WITH ALU SUPREMO MUTATE -- SOMETHING OUTSIDE OF THE MUTATE FUNCTION IS MESSING WITH THE SEQUENCE OR SOMETHIGN

plt.scatter(alu_supremo_mutate_scores['mse_mean'], random_sample_all['mse_mean'])

plt.scatter(random_sample_all['mse_mean'], random_sample_all['mse_mean_new'])
lims = [0,0, 0.032, 0.032]
plt.plot(lims, lims, 'k--', alpha=0.75, zorder=0)  # dashed black

fig, axes = plt.subplots(1, 2, figsize=(6 * 2, 5), sharex=True, sharey=True)

strands=['+', '-']
for ax, strand in zip(axes, strands):
    group = random_sample_all[random_sample_all['strand'] == strand]
    ax.scatter(group['mse_mean'], group['mse_mean_new'], alpha=0.7)
    
    lims = [0,0, 0.032, 0.032]
    ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)  # dashed black

    ax.set_title(f"Strand: {strand}")
    ax.set_xlabel("MSE mean w/ both strands as base (augmented)")
    ax.set_ylabel("MSE mean w/ + strand base only (augmented)")

plt.tight_layout()
plt.show()

random_sample_all

plt.scatter(random_sample_all['mse_mean_new'], supremo_latest_scores['mse_HFF_mean'])
lims = [0,0, 0.032, 0.032]
plt.plot(lims, lims, 'k--', alpha=0.75, zorder=0)  # dashed black

# #read in alignment sampling scores 
# alignment_sample_idx=pd.read_csv(f'{data_path}alignment_sample_indices.txt', header=None)

# final_alignment_df=pd.read_csv(f'{score_path}/combined/20240308_Alu_20alignment_sample_scores.txt', sep='\t', index_col=0)
# final_alignment_df

chrom_size_path=f'{POLLARD_DATA}/vertebrate_genomes/human/hg38/hg38/GRCh38_EBV.chrom.sizes'

chrom_sizes=pd.read_csv(chrom_size_path, sep='\t', header=None)
chrom_sizes.columns=['chrom', 'length']
chrom_sizes

Alu_sample=pd.read_csv(f'{res_dir}202520513_hg38Alus_100sample_repName_DEL_combinedscores.txt', sep='\t', index_col=0)
Alu_sample

total_length=1048576

#get coordinates of surrounding Alu window 

def get_Alu_window(row):
    
    chrom_alu=row['CHROM']
    pos=row['POS']
    posend=row['END']
    SVLEN=row['SVLEN']
    chrom_len=chrom_sizes[chrom_sizes['chrom']==chrom_alu]['length'].values[0]
    
    half=np.floor((total_length-SVLEN)/2)
    
    start = pos - half
    end = posend + half    
    
    # shift window if it goes out of bounds
    if start < 0:
        end += -start
        start = 0
    if end > chrom_len:
        shift = end - chrom_len
        start = max(0, start - shift)
        end = chrom_len

    return pd.Series({'window_start': int(start), 'window_end': int(end)})

Alu_sample[['window_start', 'window_end']] = Alu_sample.apply(get_Alu_window, axis=1)

Alu_sample

#get files with exons, genes, all Alus

repeat_dir=f'{DATA_DIR}/alus/data/repeats/'
Alus=pd.read_csv(f'{repeat_dir}hg38_repeatmaster_Alus.txt', sep='\t', index_col=0)
Alus['subFamily']=Alus['repName'].str[:4]
Alus_ALL_BED= BedTool.from_dataframe(Alus[['genoName', 'genoStart', 'genoEnd']])
Alus

gene_annot_hg38 = pd.read_csv(f'{SUPREMO_DIR}/data/gene_annot_hg38', sep = ',')
gene_annot_hg38 = gene_annot_hg38[['.' not in x for x in gene_annot_hg38.gene]]
gene_annot_hg38_BED = BedTool.from_dataframe(gene_annot_hg38[['chr', 'start', 'end', 'gene', 'strand']])
gene_annot_hg38

#these are only protein coding 
gene_annot_PC = pd.read_csv('/pollard/home/szhang20/data/hg38_gene_annot_PC', sep = ',')
gene_annot_PC = gene_annot_PC[['.' not in x for x in gene_annot_PC.gene]]
gene_annot_PC_BED = BedTool.from_dataframe(gene_annot_PC[['chr', 'start', 'end', 'gene', 'strand']])
gene_annot_PC

gtf_file=f'{POLLARD_DATA}/genetics/GENCODE_hg38/gencode.v31.annotation.gtf.gz'

# Parse the GTF file
exons = []
with gzip.open(gtf_file, 'rt') as f:
    for i, line in enumerate(f):
#         if i==20:
#             break
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        if columns[2] == 'exon':
            # Extract attributes
            info = columns[8]

            #extract gene name 
            gene_name=""
            match = re.search(r'gene_name "([^"]+)"', info)    
            if match:
                gene_name=match.group(1)  
            
            exon_chr=columns[0]
            exon_start = int(columns[3])
            exon_end = int(columns[4])
            exons.append((gene_name, exon_chr, exon_start, exon_end))

df_exons = pd.DataFrame(exons, columns=['Gene Name', 'Exon Chr', 'Exon Start', 'Exon End'])

exons_BED = BedTool.from_dataframe(df_exons[['Exon Chr', 'Exon Start', 'Exon End', 'Gene Name']])

df_exons

print(alu_BED)

alu_BED = BedTool.from_dataframe(Alu_sample[['CHROM', 'POS', 'END', 'row_idx']])
alu_window_BED= BedTool.from_dataframe(Alu_sample[['CHROM', 'window_start', 'window_end', 'row_idx']])

Alu_sample[30:40]

Alu_sample_annot=Alu_sample.copy()

feature_dict={'gene_PC' : gene_annot_PC_BED, 'gene_all' : gene_annot_hg38_BED, 
              'exons' : exons_BED, 'alus_all' : Alus_ALL_BED}
for feat, feat_BED in feature_dict.items():
    
    print(feat)

    #need to run each alu separately to calculate bp overlaps 

    counter=0
    for i, alu_region in enumerate(alu_BED):
        
        print(counter)
        counter +=1
        #print('alu region', alu_region)
        
        try:

            alu_single = BedTool([alu_region])
            intersect = alu_single.intersect(feat_BED)

            #print('intersect', intersect)
            # If no intersections, skip
            if len(intersect) == 0:
                continue

            # Sort and merge the overlaps
            merged = intersect.sort().merge()

            #print('merged', merged)

            # Compute base pair length
            nonredundant_bp_total=0
            nonredundant_bp_total += sum(f.length for f in merged)

            region_len=int(alu_region[2])-int(alu_region[1])
            #print('nonredundant_bp_total', nonredundant_bp_total)

            Alu_sample_annot.loc[i, f'{feat}_overlap']=nonredundant_bp_total/region_len
        except:
            print('alu error in ', i)

    counter=0
    for i, alu_region in enumerate(alu_window_BED):
        
        print(counter)
        counter +=1
        #print('alu region', alu_region)
        
        try:

            alu_single = BedTool([alu_region])
            intersect = alu_single.intersect(feat_BED)

            #print('intersect', intersect)
            # If no intersections, skip
            if len(intersect) == 0:
                continue

            # Sort and merge the overlaps
            merged = intersect.sort().merge()

            #print('merged', merged)

            # Compute base pair length
            nonredundant_bp_total=0
            nonredundant_bp_total += sum(f.length for f in merged)

            region_len=int(alu_region[2])-int(alu_region[1])
            #print('nonredundant_bp_total', nonredundant_bp_total)

            Alu_sample_annot.loc[i, f'{feat}_window_overlap']=nonredundant_bp_total/region_len
        except:
            print('window error in ', i)

Alu_sample_annot.to_csv('{res_dir}20250606_Alu_5100sample_overlapannot.txt', sep='\t', index=False)

Alu_alignment=pd.read_csv(f'{data_dir}20240308_hg38Alus_1000_random_alignment_sample_DEL.txt', sep='\t')
Alu_alignment=Alu_alignment[['CHROM', 'POS', 'END', 'align_percent_to_alu', 'align_percent_to_consensus']]

Alu_sample_annot=pd.merge(Alu_sample_annot, Alu_alignment, how='left')
Alu_sample_annot
Alu_sample_annot.to_csv('{res_dir}20250606_Alu_5100sample_overlapannot.txt', sep='\t', index=False)

pybedtools.cleanup()

# Create a new BedTool from the overlapping regions only
overlap_regions = A.intersect(B)

# Merge the overlapping regions to remove duplicates
merged = overlap_regions.merge()

# Sum the lengths of the merged intervals
nonredundant_bp = sum(f.length for f in merged)

# --- cell 48: did not parse; preserved as comments ---

#   pybedtools.cleanup()

# def get_exons(gene, map_start, map_end):
#     gene_exon_df_unsort=df_exons[df_exons['Gene Name']==gene]

#     gene_exon_df = gene_exon_df_unsort.sort_values(by='Exon Start')
#     gene_exon_df = gene_exon_df[gene_exon_df['Exon End'] <= map_end]
#     gene_exon_df = gene_exon_df[gene_exon_df['Exon Start'] >= map_start]

#     exon_coords=list(zip(gene_exon_df['Exon Start']-map_start, gene_exon_df['Exon End']-map_start))
#     print('exon coords', exon_coords)
#     return(exon_coords)
