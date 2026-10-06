"""Converted from alu_celltypespecific_gene.ipynb; logic unchanged."""

import sys
from collections import namedtuple
import os
import pickle
import re
import subprocess

import h5py
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import scipy.stats
import seaborn as sns
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import parse_attributes

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

#from plotnine import *

data_dir=f'{DATA_DIR}/alus/data/'

exp_table=pd.read_csv(f'{data_dir}encode/ENCFF759DYB_rnaseq_all_fromsandy.tsv', sep='\t')

#read in gencode genes to decide ones we want to filter on 
#function expand info in gencode list 

pd.set_option('display.max_columns', None)

#read in gencode information 
gencode_annot_path=f'{data_dir}gencode/gencode.v29.annotation.gtf.gz'
gencode_df = pd.read_csv(gencode_annot_path, sep='\t', comment='#', header=None,compression='gzip', low_memory=False,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "score2", "info"])

#apply to gencode info column 
parsed = gencode_df['info'].apply(parse_attributes)
gencode_df = pd.concat([gencode_df, pd.DataFrame(parsed.tolist())], axis=1)

#read in broader gene type definitions to use 
gencode_categories=pd.read_csv(f'{PROJECT_DIR}/results/gencode_df_genes_count_broad.txt', sep='\t')
gencode_categories

#df1=gencode_df_genes['gene_type'].value_counts()
#df1.to_csv(f'{PROJECT_DIR}/results/gencode_df_genes_count.txt', sep='\t')

gencode_categories

gencode_df_genes=gencode_df[gencode_df['feature']=='gene']
gencode_df_genes['gene_len']=gencode_df_genes['end']-gencode_df_genes['start']
#gencode_df_genes['gene_type'].value_counts()

#merge with the broader gene categories
gencode_df_genes=gencode_df_genes.merge(gencode_categories, on='gene_type')

#filter exp table to PC genes as in gencode 
#check if there are lncRNAs
gencode_df_PC=gencode_df_genes[gencode_df_genes['gene_type']=='protein_coding']
gencode_df_lincRNA=gencode_df_genes[gencode_df_genes['gene_type']=='lincRNA']
gencode_df_lncRNA=gencode_df_genes[gencode_df_genes['broad_gene_type']=='lncRNA']
gencode_df_pseudogene=gencode_df_genes[gencode_df_genes['gene_type'].str.contains('pseudogene')]

gencode_df_PC[gencode_df_PC['gene_name']=='TAF7']

pd.set_option('display.max_rows', None)

gencode_df[gencode_df['info'].str.contains('TAF7')]

exp_table_gene=exp_table[exp_table['gene_name'].isin(gencode_df_PC['gene_name'])]
exp_table_lincRNA=exp_table[exp_table['gene_name'].isin(gencode_df_lincRNA['gene_name'])]
exp_table_lncRNA=exp_table[exp_table['gene_name'].isin(gencode_df_lncRNA['gene_name'])]
exp_table_pseudogene=exp_table[exp_table['gene_name'].isin(gencode_df_pseudogene['gene_name'])]

def parse_cell(cell, which="left"):
    """
    parse gene exp values into tpm and fpkm values, average cross replicates 
    """
    if pd.isna(cell):
        return np.nan

    s = str(cell)

    if "_" not in s:
        # If no underscore, treat whole thing as a single value
        try:
            return float(s)
        except ValueError:
            return np.nan

    left, right = s.split("_", 1)
    part = left if which == "left" else right

    # Split by ':' and take mean of all numeric parts
    vals = []
    for x in part.split(":"):
        try:
            vals.append(float(x))
        except ValueError:
            pass

    if len(vals) == 0:
        return np.nan
    return np.mean(vals)

gene_table_dict={'protein_coding' : exp_table_gene, 'lincRNA' : exp_table_lincRNA, 'lncRNA' : exp_table_lncRNA, 'pseudogene' : exp_table_pseudogene}
gene_tpm_dict={}
gene_fpkm_dict={}

for gene_type, exp_table in gene_table_dict.items():
    print(gene_type)
    
    gene_cols = exp_table.columns[:2]
    data_cols = exp_table.columns[2:]

    # build TPM and FPKM numeric matrices
    tpm_data = exp_table[data_cols].applymap(lambda x: parse_cell(x, which="left"))
    #fpkm_data = exp_table[data_cols].applymap(lambda x: parse_cell(x, which="right"))
    
    # fix column names: keep only part before "_"
    new_cols = [c.split("_", 1)[0] for c in data_cols]
    tpm_data.columns = new_cols
    #fpkm_data.columns = new_cols
    
    # reattach gene columns
    tpm_df = pd.concat([exp_table[gene_cols], tpm_data], axis=1)
    #fpkm_df = pd.concat([exp_table[gene_cols], fpkm_data], axis=1)

    # add to dict
    gene_tpm_dict[gene_type] = tpm_df
    #gene_fpkm_dict[gene_type] = fpkm_df

mapping=pd.read_csv(f'{data_dir}encode/encode_rnaseq_experimentreport.tsv', sep='\t', skiprows=1)
acc_names=pd.DataFrame({'Accession': [i.split('_')[0] for i in exp_table.columns[2:]]})
acc_names=acc_names.merge(mapping)
#acc_names=acc_names[acc_names['Biosample classification']!='tissue']

gene_cols=['gene_id', 'gene_name']
diffcell_acc=acc_names[acc_names['Biosample classification']=='in vitro differentiated cells']['Accession'].tolist()
primarycell_acc=acc_names[acc_names['Biosample classification']=='primary cell']['Accession'].tolist()
cellline_acc=acc_names[acc_names['Biosample classification']=='cell line']['Accession'].tolist()
HFF_accession=mapping[mapping['Description'].str.contains('HFF', na=False)]['Accession'].tolist()

cell_group_dict={'diffcell' : diffcell_acc, 'primarycell' : primarycell_acc, 'cellline' : cellline_acc, 'HFF' : HFF_accession}

def row_entropy(row):
    x = row.iloc[2:].to_numpy(dtype=float) 
    x = np.nan_to_num(x, nan=0.0)
    #x = np.log2(x + 1)   
    if x.sum() == 0:
        return 0.0
    return scipy.stats.entropy(x, base=2)  

def row_simpson(row):
    """
    Calculates Simpson's Index (D).
    Higher D = Higher specificity/dominance.
    Lower D = More even/housekeeping.
    """
    x = row.iloc[2:].to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0)
    
    total = x.sum()
    if total == 0:
        return 0.0
    
    # Calculate pi (relative abundance)
    p = x / total
    # Calculate sum of squared proportions: sum(pi^2)
    return np.sum(p**2)

def row_std(row):
    x = row.iloc[2:].to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0)
    if np.all(x == 0):
        return 0.0
    return np.std(x)

def row_mean(row):
    x = row.iloc[2:].to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0)
    if np.all(x == 0):
        return 0.0
    return np.mean(x)

def row_max(row):
    x = row.iloc[2:].to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0)
    if np.all(x == 0):
        return 0.0
    return np.max(x)

gene_cols

gene_tpm_dict.keys()

# results = {}

# for gene_type, tpm_df in gene_tpm_dict.items():
#     print(gene_type)
    
#     df_out = tpm_df.copy()
    
#     for group_name, acc_cols in cell_group_dict.items():
#         expr_cols = [c for c in tpm_df.columns[2:] if c in acc_cols]
        
#         df_out[f"entropy_{group_name}"] = (df_out[expr_cols].apply(row_entropy, axis=1))
#         df_out[f"std_{group_name}"] = (df_out[expr_cols].apply(row_std, axis=1))
#         df_out[f"mean_{group_name}"] = (df_out[expr_cols].apply(row_mean, axis=1))
#         df_out[f"max_{group_name}"] = (df_out[expr_cols].apply(row_max, axis=1))

#     results[gene_type] = df_out

def apply_metrics(df, tpm_df, cell_group_dict, suffix=""):
    """Apply entropy, std, mean, max metrics to df for each cell group."""
    for group_name, acc_cols in cell_group_dict.items():
        expr_cols = [c for c in tpm_df.columns[2:] if c in acc_cols]
        col_prefix = f"{group_name}{suffix}"
        
        df[f"entropy_{col_prefix}"] = df[expr_cols].apply(row_entropy, axis=1)
        df[f"simpson_{col_prefix}"] = df[expr_cols].apply(row_simpson, axis=1)
        df[f"std_{col_prefix}"]     = df[expr_cols].apply(row_std, axis=1)
        df[f"mean_{col_prefix}"]    = df[expr_cols].apply(row_mean, axis=1)
        df[f"max_{col_prefix}"]     = df[expr_cols].apply(row_max, axis=1)
    return df

MIN_MEAN_TPM = 0   # minimum mean TPM across ALL columns
MIN_MAX_TPM  = 5   # at least one cell line must exceed this

results = {}
for gene_type, tpm_df in gene_tpm_dict.items():
    print(gene_type)

    expr_all_cols = tpm_df.columns[2:]

    # unfiltered 
    df_out = tpm_df.copy()
    df_out = apply_metrics(df_out, tpm_df, cell_group_dict, suffix="")

    # filtered
    mean_tpm = tpm_df[expr_all_cols].mean(axis=1)
    max_tpm  = tpm_df[expr_all_cols].max(axis=1)
    mask     = (mean_tpm >= MIN_MEAN_TPM) & (max_tpm >= MIN_MAX_TPM)

    df_filtered = tpm_df[mask].copy()
    df_filtered = apply_metrics(df_filtered, df_filtered, cell_group_dict, suffix="_filtered")

    # merge so both sit in one df per gene type; unmatched filtered rows get NaN in unfiltered cols
    df_out = df_out.merge(
        df_filtered[[c for c in df_filtered.columns if "_filtered" in c or c in tpm_df.columns[:2]]],
        on=list(tpm_df.columns[:2]),
        how="left"
    )

    results[gene_type] = df_out

results['protein_coding']

for gene_type, df_all in results.items():

    #filter to nonzero cols     
    cols=['max_diffcell', 'max_primarycell', 'max_cellline']
    df_out = df_all[(df_all[cols] != 0).all(axis=1)]

    ####### plot entropy 
    plt.hist(df_out['entropy_diffcell_filtered'], bins=30, alpha=0.5, label='DiffCell')
    plt.hist(df_out['entropy_primarycell_filtered'], bins=30, color='green', alpha=0.5, label='Primary')
    plt.hist(df_out['entropy_cellline_filtered'], bins=30, color='orange', alpha=0.5, label='CellLine')
    
    plt.xlabel('Entropy')
    plt.ylabel('Frequency')
    plt.title(f'{gene_type}: Filtered Entropy Across Cell Types (non-zero expression)')
    plt.legend()  
    plt.show()

    ####### plot simpson 
    plt.hist(df_out['simpson_diffcell_filtered'], bins=30, alpha=0.5, label='DiffCell')
    plt.hist(df_out['simpson_primarycell_filtered'], bins=30, color='green', alpha=0.5, label='Primary')
    plt.hist(df_out['simpson_cellline_filtered'], bins=30, color='orange', alpha=0.5, label='CellLine')
    
    plt.xlabel('Entropy')
    plt.ylabel('Frequency')
    plt.title(f'{gene_type}: Filtered Simpson Across Cell Types (non-zero expression)')
    plt.legend()  
    plt.show()

    ####### compare entropy vs simpson 
    plt.scatter(df_out['entropy_diffcell_filtered'], df_out['simpson_diffcell_filtered'], alpha=0.5, label='DiffCell')
    plt.scatter(df_out['entropy_primarycell_filtered'], df_out['simpson_primarycell_filtered'], color='green', alpha=0.5, label='Primary')
    plt.scatter(df_out['entropy_cellline_filtered'], df_out['simpson_cellline_filtered'], color='orange', alpha=0.5, label='CellLine')
    
    plt.xlabel('Entropy')
    plt.ylabel('Simpson')
    plt.title(f'{gene_type}')
    plt.legend()  
    plt.show()

# plt.hist(tpm_df_diffcell['std_diffcell'], bins=1000, alpha=0.5)
# plt.hist(tpm_df_primary['std_primary'], bins=1000, color='green', alpha=0.5)
# plt.hist(tpm_df_cellline['std_cellline'], bins=1000, color='orange', alpha=0.5)
# plt.xlim(0,100)

#after pulling out entropy
#plot entropy of gene vs #Alus / gene len, # high disruptive Alus / gene len, entropy of Alu scores within gene 
#better to use SD since not normal distributions

gencode_df_genes

#os.listdir(alu_res_dir)

# alu_res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
# alu_annot_withtracks_het=pd.read_csv(f'{alu_res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', 
#                                      sep='\t', index_col=0)

alu_res_dir=f'{PROJECT_DIR}/results/paper_results/'
alu_annot_withtracks_het=pd.read_csv(f'{alu_res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt', 
                                     sep='\t', index_col=0)

alu_annot_withtracks_het['single_gene_name'].value_counts()[0:20]

(1109106-422114)/1109106

#merge back gene ids just in case, so can merge them with the encode rna seq
alu_annot_withtracks_het=alu_annot_withtracks_het.merge(gencode_df_genes[['gene_id', 'gene_name']], left_on='single_gene_name', right_on='gene_name')

alu_annot_withtracks_het

no_na_alu=alu_annot_withtracks_het.dropna(subset='mse_mean')
p99 = no_na_alu["mse_mean"].quantile(0.99)
p50 = no_na_alu["mse_mean"].quantile(0.50)

df_top1 = no_na_alu[no_na_alu["mse_mean"] >= p99]

df_top1

1-(3484/11092)

df_top1['single_gene_name'].value_counts()

agg_alu_gene={}

#merge agg_df with entropy information 
for gene_type, df_out in results.items():
    print(gene_type)

    #filter out df to ones that overlap the genes we want 
    df_overlap = alu_annot_withtracks_het[(alu_annot_withtracks_het["single_gene_name"]!='None')]
    df_overlap=df_overlap.merge(gencode_categories, how='left', left_on='single_gene_type', right_on='gene_type')

    if gene_type=='lncRNA':
        df_valid=df_overlap[df_overlap['broad_gene_type']=='lncRNA']
    elif gene_type=='pseudogene':
        df_valid=df_overlap[df_overlap['single_gene_type'].str.contains('pseudogene')]
    else:
        df_valid=df_overlap[df_overlap['single_gene_type']==gene_type]

    print(len(df_valid))
    # total number of rows per gene
    total_counts = df_valid.groupby("single_gene_name")["single_gene_name"].count()
    #high MSE rows
    high_MSE = df_valid[df_valid["mse_mean"] > p99].groupby("single_gene_name")["single_gene_name"].count()
    #average MSE vals
    avg_mse = df_valid.groupby("single_gene_name")["mse_mean"].mean()
    std_mse = df_valid.groupby("single_gene_name")["mse_mean"].std()

    agg_df = pd.DataFrame({
        "gene": total_counts.index,
        "total_alu": total_counts.values,
        "high_alu": high_MSE.reindex(total_counts.index).fillna(0).values,
        "avg_alu_mse": avg_mse.reindex(total_counts.index).values, 
        "std_alu_mse": std_mse.reindex(total_counts.index).values  
    })
    
    agg_df["high_alu_gene_dens"] = agg_df["high_alu"] / agg_df["total_alu"]
    
    # join gene lengths
    agg_df = agg_df.merge(gencode_df_genes[["gene_name", "gene_id", "gene_len", 'gene_type']],
        left_on="gene", right_on="gene_name", how="left")
        #on="gene_id", how="left")
    
    # alu density per gene length
    agg_df["alu_dens_per_len"] = 100*(agg_df["total_alu"] / agg_df["gene_len"])
    agg_df["alu_avg_std_divided_mean"] = (agg_df["std_alu_mse"] / agg_df["avg_alu_mse"]) 
    
    agg_df = agg_df[["gene", "gene_id", 'gene_type', "high_alu_gene_dens", "alu_dens_per_len", 'high_alu', 'total_alu', 
                     'avg_alu_mse', 'std_alu_mse','alu_avg_std_divided_mean', 'gene_len']]

    #finally, merge with entropy information 
    cols_want=[]
    prefixes = ('gene_name', 'entropy', 'simpson', 'std', 'mean', 'max')  
    cols_want = df_out.columns[df_out.columns.str.startswith(prefixes)]
    agg_df_expentropy=agg_df.merge(df_out[cols_want], left_on='gene', right_on='gene_name', how='left')

    agg_alu_gene[gene_type]=agg_df_expentropy

agg_alu_gene['protein_coding']

(agg_alu_gene['protein_coding']['entropy_cellline_filtered']).isna().sum()

#save each into a df
for gene_type, df_agg in agg_alu_gene.items():  
    df_agg.to_csv(f'{PROJECT_DIR}/results/20260218_CTS_genes/20260406_{gene_type}_entropysimpson_TPM5.txt', sep='\t')

    #need to also merge back to original Alu table so we know which Alu transcribed

agg_alu_gene

for gene_type, df_agg in agg_alu_gene.items():

    #filter to nonzero cols      
    cols=['entropy_cellline_filtered', 'entropy_primarycell_filtered', 'entropy_diffcell_filtered']
    cols=['simpson_cellline_filtered', 'simpson_primarycell_filtered', 'simpson_diffcell_filtered']

    for col in cols:
        #plot # of high MSE Alus vs entropy 
        plt.scatter(df_agg['avg_alu_mse'], df_agg[col])

        #plot std of scores vs entropy 
        #plt.scatter(df_agg[col], df_agg['avg_alu_mse'])
        
        plt.xlabel('# of High MSE Alus in Gene')       # x-axis label
        plt.ylabel(f'{col}') 
        plt.title(f'{gene_type}: {col} vs # of High MSE Alus')
        #plt.xscale('log')
        plt.show()

for gene_type, df_agg in agg_alu_gene.items():

    cols = ['entropy_cellline', 'entropy_primarycell', 'entropy_diffcell']

    for col in cols:

        df = df_agg[['avg_alu_mse', col]].dropna().copy()

        # bin x-axis into quantiles
        df['mse_bin'] = pd.qcut(df['avg_alu_mse'], q=10, duplicates='drop')

        # prepare data for boxplot
        grouped = [g[col].values for _, g in df.groupby('mse_bin')]

        plt.boxplot(grouped, showfliers=False)

        plt.xticks(
            range(1, len(df['mse_bin'].unique()) + 1),
            [str(b) for b in df['mse_bin'].unique()],
            rotation=45
        )

        plt.xlabel('avg_alu_mse quantile bins')
        plt.ylabel(col)
        plt.title(f'{gene_type}: {col} vs avg_alu_mse quantiles')

        plt.tight_layout()
        plt.show()

high_alu_list=[]
for gene_type, df_agg in agg_alu_gene.items():
    high_alu=df_agg[df_agg['high_alu']>20]
    high_alu_list.append(high_alu)
high_alu_df=pd.concat(high_alu_list, axis=0)

len(high_alu_df['gene'].tolist())

high_alu_df

non_string_cols = agg_df_expentropy.select_dtypes(exclude="object").columns.tolist()

corr = agg_df_expentropy[non_string_cols].corr(method='pearson')

plt.figure(figsize=(10, 8))
sns.heatmap(corr, annot=True, fmt=".2f", cmap='coolwarm', center=0)
plt.title("Genes with Alus vs Expression Entropy")
#plt.savefig(f"{PROJECT_DIR}/Alu_feature_overlap_1kb_up.png", dpi=300, bbox_inches="tight", transparent=True)
plt.show()

agg_df_expentropy

plt.scatter(agg_df_expentropy['alu_dens_per_len'], agg_df_expentropy['entropy_diffcell'])
plt.xlabel('High # Alu')       # x-axis label
plt.ylabel('Entropy (Primary)')

plt.scatter(agg_df_expentropy['avg_alu_mse'], agg_df_expentropy['entropy_diffcell'])
plt.xlabel('Average Alu DEL MSE in gene')       # x-axis label
plt.ylabel('Entropy (Diff Cell)') 
#plt.xscale('log')

# ## Enformer targets
# #### explore heterochromatin

#parse enformer targets
#human_targets=pd.read_csv(f'{data_dir}enformer_targets_human.txt', sep='\t')

#use new one with modified description 
human_targets_file = f'/pollard/home/szhang20/SuPreMo/Enformer_model/targets_human_descriptionindex.txt'
human_targets=pd.read_csv(human_targets_file, sep='\t')
human_targets['assay']=human_targets['description'].str.split(':').str[0]
human_targets[human_targets['description'].str.contains('H3K27me3')]

#to get encode tsv information
#under 'Experiment Search', need to use 'List' as sort, and then download tsv 
chip_exp_metadata=pd.read_csv(f'{data_dir}encode/encode_chip_experimentreport.tsv', sep='\t', skiprows=1)

#combine with ENCODE metadata information (string is under 'Files') 

# Build one regex that matches any human_targets value
identifiers = human_targets['identifier'].astype(str).tolist()
pattern = '(' + '|'.join(map(re.escape, identifiers)) + ')'

# Extract the matching target from 'files'
chip_exp_metadata['identifier'] = chip_exp_metadata['Files'].str.extract(pattern, expand=False)

#merge
cols_tomerge=chip_exp_metadata.columns[0:20].tolist() + ['identifier']
merged = human_targets.merge(chip_exp_metadata[cols_tomerge], on='identifier', how='left')

chip_enformer=merged[merged["ID"].notna()]
chip_enformer

#marks i want

#constitutive heterochromatin 
#H3K9me3
#TRIM28 (recruited by KRABZNF)
#SETDB1 (recruited by KRABZNF)
#CBX proteins (part of HP1), binds to H3K9me3
#DNMT1

#Facultatitve heterochromatin 
#H3K9me2, H3K27me3
#EZH2, SUZ12, RNF2 #components of polycomb complex 
#REST, YY1

#BDP1

#anything that starts with ZNF

list_to_contain=['H3K9me3', 'H3K9me2', 'H3K27me3', 'TRIM28', 'SETDB1', 'CBX', 'DNMT1', 'EZH2', 'SUZ12', 'RNF2', 'REST', 'YY1', 'ZNF', 'CTCF']
# Build a regex pattern that matches any of the substrings
pattern_tocontain = '|'.join(list_to_contain)

#filter to those with replicates, and cell lines/primary cells  
filtered = chip_enformer[
    chip_enformer['Target of assay'].str.contains(pattern_tocontain) &
    chip_enformer['Biological replicate'].str.contains(',') &
    (chip_enformer['Biosample classification'] != 'tissue') & 
    (~chip_enformer['description'].str.lower().str.contains('genetically modified'))
]

all_targets_str = ' '.join(f'"{str(x)}"' for x in filtered['description'])
all_targets_str

#filtered[filtered['description'].str.contains('CTCF')]['description'].tolist()

alu_res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'

#finally, let's pick the ones we want to test 
alu_res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
alu_annot_withtracks_het=pd.read_csv(f'{alu_res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', 
                                     sep='\t', index_col=0)

alu_annot_withtracks_het

test1=alu_annot_withtracks_het.sample(20)

#split out cols to come first 
col_first=['CHROM', 'POS', 'REF', 'ALT', 'END', 'SVTYPE', 'SVLEN', 'strand', 'orig_idx']
remaining_columns = [col for col in test1.columns if col not in col_first]
col_order = col_first + remaining_columns

test1=test1[col_order]
test1

test2=alu_annot_withtracks_het[['CHROM', 'REF_start', 'REF_stop']]
#test2['REF_start']=test2['REF_start'].astype(int)
#test2['REF_stop']=test2['REF_stop'].astype(int)

df = test2.copy().dropna(subset=['REF_start','REF_stop'])
df['REF_start']=df['REF_start'].astype(int)
df['REF_stop']=df['REF_stop'].astype(int)

df.copy().dropna(subset=['REF_start','REF_stop'])

df[df['REF_stop']=='']

test2=alu_annot_withtracks_het[['CHROM', 'REF_start', 'REF_stop']]
test2['REF_start']=test2['REF_start'].astype(int)
test2['REF_stop']=test2['REF_stop'].astype(int)

#separate again into top 1 and bottom 50%

test1.to_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/alu_test.txt', sep='\t')

test2=pd.read_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/alut_test_scores', sep='\t')
test2.columns

selected_tracks=["CHIP:H3K27me3:HepG2", "CHIP:H3K27me3:HeLa-S3", "CHIP:H3K9me3:K562"]

human_targets_file = f'/pollard/home/szhang20/SuPreMo/Enformer_model/targets_human_descriptionindex.txt'

human_targets = pd.read_csv(human_targets_file, sep='\t')
#human_targets["description"] = human_targets["description"] + ":" + human_targets["index"].astype(str)
#human_targets.to_csv('/pollard/home/szhang20/SuPreMo/Enformer_model/targets_human_descriptionindex.txt', sep='\t', index=False)

#read in annotated/scored Alus 

#let's say we want to shuffle 5000 of the top 1% 
#the nanother 5000 randomly selected from the rest of the genome
#for each alu we shuffle 10x with regions nerably --> total 100,000 shuffles 

alu_scores=alu_annot_withtracks_het.dropna(subset='mse_mean')

mse_1perc= alu_scores['mse_mean'].quantile(0.99)
top1_mse = alu_scores[alu_scores['mse_mean'] >= mse_1perc]
top1_mse_sample=top1_mse.sample(n=5000, replace=False, random_state=729)

mse_50perc= alu_scores['mse_mean'].quantile(0.50)
bottom50_mse=alu_scores[alu_scores['mse_mean'] <= mse_50perc]
bottom50_mse_sample=bottom50_mse.sample(n=5000, replace=False, random_state=729)

#alu_to_shuffle=pd.concat([top1_mse_sample, bottom50_mse_sample], axis=0)

top1_test=top1_mse_sample.sample(n=10, replace=False, random_state=729)
top1_test.to_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/test/top1mse_50test.txt', sep='\t')

# top1_mse_sample.to_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/top1mse_5000samp.txt', sep='\t')
# bottom50_mse_sample.to_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/bottom50mse_5000samp.txt', sep='\t')

alu_to_shuffle

#ok, let's read in the scores 
top1_scores=pd.read_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/top1mse_5000samp_scores', sep='\t')
bottom50_scores=pd.read_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/bottom50mse_5000samp_scores', sep='\t')

top1_scores

mse_cols = top1_scores.filter(regex=r'^corr.*0$').columns
top1_scores_subset = top1_scores[mse_cols].reset_index()
top1_mse_annot=top1_mse_sample[['CHROM', 'orig_idx', 'corr_mean', 'HFF_H3K9me3_window_overlap', 'HFF_H3K27me3_window_overlap']].reset_index()
top_all=pd.concat([top1_mse_annot, top1_scores_subset], axis=1)
top_all

bottom50_scores_subset = bottom50_scores[mse_cols].reset_index()
bottom50_annot=bottom50_mse_sample[['CHROM', 'orig_idx', 'corr_mean', 'HFF_H3K9me3_window_overlap', 'HFF_H3K27me3_window_overlap']].reset_index()
bottom_all=pd.concat([bottom50_annot, bottom50_scores_subset], axis=1)

test2=pd.concat([top_all, bottom_all], axis=0)
test2

#top1_mse_sample.sort_values(by='mse_mean', ascending=False)[0:50]

corrs = test2.corr()['corr_mean']
corrs.sort_values(ascending=False)[0:50]

corrs.loc[corrs.index.str.contains('RNF2')]

#loop through and plot 
corrs = test2.corr()['mse_mean'].drop('mse_mean', errors='ignore')
top25 = corrs.sort_values(ascending=False).head(25)

for col in top25.index:
    plt.figure(figsize=(5, 4))
    
    plt.hist(bottom_all[col],bins=30, color='blue', alpha=0.5, label='bottom50')
    plt.hist(top_all[col],bins=30, color='red', alpha=0.5, label='top1')

    plt.xlabel(col)
    plt.ylabel('Count')
    plt.legend()
    
    plt.tight_layout()
    plt.show()

plt.hist(bottom_all['mse_CHIP:TRIM28:K562:1053_1053_0'], bins=30, color='blue', alpha=0.5)
plt.hist(top_all['mse_CHIP:TRIM28:K562:1053_1053_0'], bins=30, color='red', alpha=0.5)

plt.scatter( top_all['mse_CHIP:ZNF579:MCF-7:1504_1504_0'], top_all['mse_mean'], c=top_all['HFF_H3K27me3_window_overlap'], s=20, alpha=0.9)
plt.yscale('log')

plt.hist(bottom_all['mse_CHIP:ZNF830:K562:4079_4079_0'],alpha=0.5)
plt.hist(top_all['mse_CHIP:ZNF830:K562:4079_4079_0'], color='red', alpha=0.5)

corrs = top_all.corr()['corr_mean']

corrs.sort_values(ascending=False)[0:20]

corrs.sort_values(ascending=True)[0:20]

corr_cols = top1_scores.filter(regex=r'^corr.*0$').columns
top1_scores_corrsubset = top1_scores[corr_cols]
test2=pd.concat([top1_mse_sample[['orig_idx', 'corr_mean', 'HFF_H3K9me3_window_overlap', 'HFF_H3K27me3_window_overlap']], 
                 top1_scores_corrsubset], axis=1)

corrs = test2.corr()['corr_mean']

corrs.sort_values(ascending=False)[0:20]

#let's look at some specific tracks

file_path = f'{PROJECT_DIR}/results/20260205_enformer_het/test/top1mse_50test_tracks.pkl'

# Open the file in binary read mode ('rb')
with open(file_path, 'rb') as file:
    data = pickle.load(file)

# Now 'data' contains the Python object (e.g., list, dictionary, model)
print(data)

variant_tracks_all = {}
with h5py.File(f'{PROJECT_DIR}/results/20260205_enformer_het/test/top1mse_50test_tracks.h5', 'r') as f:
    for key in f.keys():
        variant_tracks_all[key] = np.array(f[key])  # converts h5 dataset to ndarray

test_var=pd.read_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/test/top1mse_50test.txt', sep='\t')
test_scores=pd.read_csv(f'{PROJECT_DIR}/results/20260205_enformer_het/test/top1mse_50test_scores', sep='\t')

test_all=pd.concat([test_var, test_scores], axis=1)
test_all

test_scores

# --- cell 117: did not parse; preserved as comments ---

# mse_mean                                      1.000000
# mse_CHIP:RNF2:K562:1883_1883_0                0.535306
# mse_CHIP:DNMT1:K562:4476_4476_0               0.507629
# mse_CHIP:SETDB1:K562:1394_1394_0              0.504296
# mse_CHIP:CBX3:K562:1052_1052_0                0.500344
# mse_CHIP:ZNF579:MCF-7:1504_1504_0             0.492790
# mse_CHIP:ZNF184:K562:3368_3368_0              0.488133
# mse_CHIP:ZNF639:K562:2989_2989_0              0.469256
# mse_CHIP:ZNF184:K562:3125_3125_0              0.468463
# mse_CHIP:ZNF207:HepG2:3558_3558_0             0.459957
# mse_CHIP:CBX1:K562:4357_4357_0                0.457527
# mse_CHIP:ZNF687:MCF-7:4223_4223_0             0.450972
# mse_CHIP:RNF2:HepG2:3077_3077_0               0.447370
# mse_CHIP:ZNF318:K562:2527_2527_0              0.446181
# mse_CHIP:ZNF263:HEK293:1381_1381_0            0.441861
# mse_CHIP:ZNF207:GM12878:1824_1824_0           0.432354
# mse_CHIP:ZNF217:GM12878:3807_3807_0           0.421852
# mse_CHIP:ZNF24:HepG2:3896_3896_0              0.419089
# mse_CHIP:RNF2:K562:3331_3331_0                0.418814
# mse_CHIP:ZNF282:K562:3746_3746_0              0.416563
# mse_CHIP:ZNF687:GM12878:4102_4102_0           0.413770
# mse_CHIP:ZNF592:MCF-7:3612_3612_0             0.413567
# HFF_H3K27me3_window_overlap                   0.413444
# mse_CHIP:ZNF384:HepG2:1767_1767_0             0.411984
# mse_CHIP:ZNF512B:MCF-7:3800_3800_0            0.411971
# mse_CHIP:ZNF384:GM12878:1246_1246_0           0.390883

col_to_test=test_scores.columns[(test_scores.columns.str.contains('DNMT1')) & (test_scores.columns.str.contains('mse'))]
test_scores[col_to_test]

df=test_scores
mse_cols = [col for col in df.columns if col.startswith('mse')]
mse_df = df[mse_cols]

#get 
flat_values = mse_df.values.flatten()
flat_indices = np.argsort(flat_values)[::-1]  

top_n = 10
top_indices = flat_indices[:top_n]

results = []
n_cols = mse_df.shape[1]

for idx in top_indices:
    row_idx, col_idx_in_mse = divmod(idx, n_cols)
    col_name = mse_df.columns[col_idx_in_mse]
    val = mse_df.values[row_idx, col_idx_in_mse]
    results.append((row_idx, col_name, val))

# 4. Print results
for r in results:
    print(f"row: {r[0]}, col: {r[1]}, mse: {r[2]}")

matches = [k for k in variant_tracks_all if 'CHIP:H3K9me3:AG04450' in k]
print(matches)

def plot_tracks(tracks, interval, height=3):
    fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True, sharey=False)
    for ax, (title, y) in zip(axes, tracks.items()):
        ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)

        #horizontal lines
        max_val = max(y)
        n_lines=5
        line_vals = np.linspace(0, max_val, n_lines + 1)[1:]  # skip 0
        for val in line_vals:
            ax.axhline(val, color='gray', linestyle=':', linewidth=0.7)

    ax.set_xlabel(f"{interval.start}-{interval.end}")
    plt.tight_layout()
    plt.show()

Interval = namedtuple("Interval", ["start", "end"])
bin_size=128
midpoint = 196450 * bin_size
interval = Interval(midpoint - 28672, midpoint + 28672)  # 448 bins * 128bp

REF_track=variant_tracks_all['1_track_CHIP:H3K9me3:AG04450:1105_REF_vrp196607_0']
ALT_track=variant_tracks_all['1_track_CHIP:H3K9me3:AG04450:1105_ALT_vrp196456_0']
DIFF_track=ALT_track - REF_track 
tracks = {'REF': REF_track, 
        'ALT': ALT_track, 
        'DIFF': DIFF_track 
          
         }
plot_tracks(tracks, interval)

#look at windows with CCRE
#go from akita side and start with high disruption and in heterochromatin alus

Interval = namedtuple("Interval", ["start", "end"])
bin_size=128
midpoint = 196450 * bin_size
interval = Interval(midpoint - 28672, midpoint + 28672)  # 448 bins * 128bp

REF_track=variant_tracks_all['9_track_CHIP:YY1:K562:973_REF_vrp196607_0_revcomp']
ALT_track=variant_tracks_all['9_track_CHIP:YY1:K562:973_ALT_vrp196459_0_revcomp']
DIFF_track=ALT_track - REF_track 
tracks = {'REF': REF_track, 
        'ALT': ALT_track, 
        'DIFF': DIFF_track 
          
         }
plot_tracks(tracks, interval)
