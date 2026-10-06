"""Converted from 20250206_alu_feature_elasticnet.ipynb; logic unchanged."""

from functools import reduce
import os
import re, os 
import subprocess

from pybedtools import BedTool
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.inspection import permutation_importance
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.metrics import roc_auc_score, average_precision_score, classification_report
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, classification_report
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import pyBigWig
import seaborn as sns
import shap

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

os.environ["OMP_NUM_THREADS"] = "12"
os.environ["OPENBLAS_NUM_THREADS"] = "12"
os.environ["MKL_NUM_THREADS"] = "12"
os.environ["VECLIB_MAXIMUM_THREADS"] = "12"
os.environ["NUMEXPR_NUM_THREADS"] = "12"

"""
###### proportion overlap 
#repeatcount 
repeatmasker=repeat_all[['genoName','genoStart', 'genoEnd']]
#gene_PC
gencode_df_genes[gencode_df_genes['gene_type']=='protein_coding'][['chrom', 'start', 'end', 'gene_name']]
#pseudogene
pseudogene=gencode_df_genes[gencode_df_genes['gene_type'].str.contains('pseudogene')][['chrom', 'start', 'end', 'gene_name']]
#lincRNA
lincRNA=gencode_df_genes[gencode_df_genes['gene_type'].str.contains('lincRNA')][['chrom', 'start', 'end', 'gene_name']]
#ZNF
ZNF=gencode_df_genes[gencode_df_genes['gene_name'].str.contains('ZNF')][['chrom', 'start', 'end', 'gene_name']]

#H3K9me3
H3K9me3=pd.read_csv(f'{data_dir}encode/HFF_H3K9me3_ENCFF327WJS.bed.gz', sep='\t', header=None)
#H3K27me3
H3K27me3=pd.read_csv(f'{data_dir}encode/HFF_H3K27me3_ENCFF979AVA.bed.gz', sep='\t', header=None)
#CTCF
HFF_CTCF=pd.read_csv(f'{data_dir}encode/HFF_CTCF_ENCFF294RSZ.bed.gz', sep='\t', header=None)

###### total count 
#alu count
#promoter -- think of another way to represent this data? maybe just counts? 
promoter=promoter_df[['chrom', 'start', 'end']]
#WGBS

######alu specific
#alignment to consensus
#consensus alignment to alu 
#divergence scores (from repeatmasker)

H3K9me3=pd.read_csv(f'{data_dir}encode/HFF_H3K9me3_ENCFF327WJS.bed.gz', sep='\t', header=None)
H3K27me3=pd.read_csv(f'{data_dir}encode/HFF_H3K27me3_ENCFF979AVA.bed.gz', sep='\t', header=None)
WGBS=pd.read_csv(f'{data_dir}encode/H1_WGBS_ENCFF434CNG.bed.gz', sep='\t', header=None)

"""

OUT_DIR = f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'
RES_DIR = f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
DATA_DIR = f'{DATA_DIR}/alus/data/'
os.listdir(OUT_DIR)

alu_annot = pd.read_csv(f'{RES_DIR}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)

#remove alus where the window contains

alu_annot.iloc[:,0:15]

alu_annot['CHROM'].value_counts()

encode_blacklist=pd.read_csv(f'{DATA_DIR}/alus/data/encBlacklist.bed', 
                             sep='\t', header=None, names=['chrom', 'start', 'end', 'reason'])
encode_blacklist['region_len']=encode_blacklist['end']-encode_blacklist['start']
encode_blacklist

encode_blacklist[0:50]

alu_annot['REF_start']=alu_annot['REF_start'].astype(int)
alu_annot['REF_stop']=alu_annot['REF_stop'].astype(int)

alu_annot_BED=BedTool.from_dataframe(alu_annot[['CHROM', 'REF_start', 'REF_stop', 'orig_idx']])
encblacklist_BED=BedTool.from_dataframe(encode_blacklist[['chrom', 'start', 'end']])

# f=0.05 means 5% of A must overlap B
# v=True removes those overlaps (invert match)
result = alu_annot_BED.intersect(encblacklist_BED, v=True, f=0.05)

# Save or use the result
alu_noblacklist_df=result.to_dataframe()

alu_annot_noblacklist_df=alu_annot[alu_annot['orig_idx'].isin(alu_noblacklist_df['name'])]
alu_annot_noblacklist_df

#do another overlap for if Alu is within a blacklist region 

alu_annot_noblacklist_df['POS']=alu_annot_noblacklist_df['POS'].astype(int)
alu_annot_noblacklist_df['END']=alu_annot_noblacklist_df['END'].astype(int)

alu_annot_noblacklist_BED=BedTool.from_dataframe(alu_annot_noblacklist_df[['CHROM', 'POS', 'END', 'orig_idx']])
encblacklist_BED=BedTool.from_dataframe(encode_blacklist[['chrom', 'start', 'end']])

# f=0.05 means 5% of A must overlap B
# v=True removes those overlaps (invert match)
result2 = alu_annot_noblacklist_BED.intersect(encblacklist_BED, v=True, f=0.05)

# Save or use the result
alu_noblacklist_df2=result2.to_dataframe()

alu_annot_noblacklist_df2=alu_annot[(alu_annot['orig_idx'].isin(alu_noblacklist_df2['name']))]
alu_annot_noblacklist_df2

OUT_DIR

# ### Add Annotations!

os.listdir(OUT_DIR)

#do a check for phyloP
test1=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_short_-2_2.csv')
test2=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_short_0_2.csv')

phyloP100=f'{POLLARD_DATA}/wynton/consortia/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw'
bw = pyBigWig.open(phyloP100)

vals = bw.values('chr1', int(7864179.0), int(8912755.0), numpy=True)  # return numpy array directly

#H1 ones -- need to change names 
CTCF_H1=pd.read_csv(f'{OUT_DIR}aluhg38_feature_count_bed_H1_CTCF.csv')
feat_overlap_H1=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_H1_histone.csv')
methylation=pd.read_csv(f'{OUT_DIR}aluhg38_WGBS_count.csv')

#replace names in H1 columns
CTCF_H1.columns = CTCF_H1.columns.str.replace('HFF', 'H1')
feat_overlap_H1.columns = feat_overlap_H1.columns.str.replace('HFF', 'H1')

#keep only IMR90 methylation column 
methylation_IMR90 = methylation.loc[:, methylation.columns.str.startswith(('CHROM', 'POS', 'END', 'orig_idx', 'IMR90'))]

methylation_IMR90

bw_avg=pd.read_csv(f'{OUT_DIR}aluhg38_feature_averageBW.csv')
GC_avg=pd.read_csv(f'{OUT_DIR}aluhg38_feature_GC.csv')
feat_overlap=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_bed.csv')
feat_count=pd.read_csv(f'{OUT_DIR}aluhg38_feature_count_bed.csv')
repeat_1mb=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_repeatsonly.csv')
repeat_other=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_repeatsonly_short.csv')

#NEW ONES
CTCF_shared=pd.read_csv(f'{OUT_DIR}aluhg38_feature_count_CTCF_sharedpeaks.csv')
segdup_selfchain=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_segdup_selfchain.csv')
phyloP_1Mb=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_1Mb_-2_2.csv')
phyloP_other=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_short_-2_2.csv')

#merge original alu with all annotations 
alu_annot_subset=alu_annot[['CHROM', 'POS', 'END', 'repName', 'strand', 'orig_idx', 'mse_mean', 'corr_mean', 
                            #'REF_start', 'REF_stop', 'var_rel_pos_REF', 
                            'single_gene_name', 'single_gene_type', 'single_gene_overlap_bp',
                           'swScore', 'milliDel', 'id','#bin', 'milliIns', 'milliDiv', 'SVLEN']]
df_list=[alu_annot_subset, GC_avg, bw_avg, feat_count, feat_overlap, repeat_other, repeat_1mb, 
         phyloP_other, phyloP_1Mb, segdup_selfchain, CTCF_H1, feat_overlap_H1, methylation_IMR90
         #CTCF_shared
        ]
merge_keys=['CHROM', 'POS', 'END','orig_idx']
def safe_merge(l, r):
    # drop any columns in r that already exist in l (except keys)
    duplicate_cols = [col for col in r.columns if col in l.columns and col not in merge_keys]
    if duplicate_cols:
        print(f"Dropping duplicate cols: {duplicate_cols}")  # helpful for debugging
        r = r.drop(columns=duplicate_cols)
    return pd.merge(l, r, on=merge_keys, how='outer')

merged = reduce(safe_merge, df_list)
#merged = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END','orig_idx'], how="outer"), df_list)

merged

# test1=merged[merged['mse_mean']>0.0286]
# sc=plt.scatter(merged['mapp36_1Mb'], merged['mapp100_1Mb'], c=np.log10(merged['mse_mean']), s=5)
# plt.colorbar(sc, label='log10(MSE mean)')
# plt.xlabel('mapp36_1Mb')
# plt.ylabel('mapp100_1Mb')

#merged.to_csv(f'{PROJECT_DIR}/results/paper_results/20260311_merged_aluhg38_all_featureannot.txt', sep='\t')

# alu_filt=merged[merged['orig_idx'].isin(alu_noblacklist_df['name'])]
# no_na_alu=alu_filt.dropna(subset='mse_mean')
# no_na_alu
#no_na_alu.to_csv(f'{PROJECT_DIR}/results/paper_results/20260311_merged_aluhg38_all_featureannot_NA_blacklist_filter.txt', sep='\t')

alu_filt_blacklist_alu_region=merged[merged['orig_idx'].isin(alu_noblacklist_df2['name'])]
no_na_alu=alu_filt_blacklist_alu_region.dropna(subset='mse_mean')
no_na_alu

#no_na_alu.to_csv(f'{PROJECT_DIR}/results/paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt', sep='\t')
#no_na_alu.to_csv(f'{PROJECT_DIR}/results/paper_results/20260316_merged_aluhg38_all_featureannot_withCTCF_NA_blacklist_filter_region_alu.txt', sep='\t')
no_na_alu.to_csv(f'{PROJECT_DIR}/results/paper_results/20260316_merged_aluhg38_all_featureannot_withH1_NA_blacklist_filter_region_alu.txt', sep='\t')

# ## Linear Regression

df_regress=no_na_alu

#filter out columns that contain phastCon or phyloP30, and also filter out phyloP bigwig columns
#also filter out the merged CTCF peaks since those are for later 
substrings_common = ['mse_mean', 'CHROM', 'orig_idx', 'GC_', 'mapp36', 'mapp100', 'phastCon100', 'alu_count', 'promoter_count', 
                     'gene_PC_overlap', 'pseudogene_overlap', 'lincRNA_overlap', 'ZNF_overlap', 'repeats_overlap',
                     'phyloP100_highfrac', 'phyloP100_lowfrac', 'seg_dup_overlap', 'self_chain_overlap',
                    #'constitutive', 'broadly_shared', 'lineage_restricted', 'cell_type_specific', 'singleton'
                     ]

substrings_HFF=['HFF_CTCF', 'HFF_H3K9me3', 'HFF_H3K27me3', 'IMR90_']
substrings_H1=['H1_CTCF', 'H1_H3K9me3', 'H1_H3K27me3', 'H1ESC_WGBS']

pattern_common = '|'.join(substrings_common)
pattern_HFF    = '|'.join(substrings_HFF)
pattern_H1     = '|'.join(substrings_H1)

common_mask = df_regress.columns.str.contains(pattern_common, case=False)
hff_mask = df_regress.columns.str.contains(pattern_HFF, case=False)
h1_mask  = df_regress.columns.str.contains(pattern_H1,  case=False)

df_HFF = df_regress.loc[:, common_mask | hff_mask]
df_H1  = df_regress.loc[:, common_mask | h1_mask]

df_HFF.columns.values

cols_1Mb = df_regress.columns[df_regress.columns.str.endswith("1Mb")]
cols_100kb = df_regress.columns[df_regress.columns.str.endswith("100kb")]
cols_1kb = df_regress.columns[df_regress.columns.str.endswith("1kb")]
cols_alu = df_regress.columns[df_regress.columns.str.endswith("alu")]

#feature_list=
feature_list=cols_1Mb.tolist()
print(feature_list)

EN_dir=f'{PROJECT_DIR}/results/20260212_elasticnet/'
CT_dict={'HFF' : df_HFF, 'H1': df_H1}

for cell_type in CT_dict:
    
    df_regress=CT_dict[cell_type]

    final_df=[]
    for feat_range in ['alu', '1kb', '10kb', '100kb', '1Mb']:
    #for feat_range in ['1Mb']:
        print(f'running {feat_range}...')

        cols_feat_range = df_regress.columns[df_regress.columns.str.endswith(feat_range)]
        feature_list=cols_feat_range.tolist()
        #feature_list=feature_list_subset

        #replace nas in feature list with 0
        df_regress_nona=df_regress.dropna(subset='mse_mean')
        df_regress_nona[feature_list] = df_regress_nona[feature_list].fillna(0)

        test_chr='chr2'

        train_df = df_regress_nona[df_regress_nona['CHROM'] != test_chr].copy()
        test_df  = df_regress_nona[df_regress_nona['CHROM'] == test_chr].copy()

        X_train = train_df[feature_list].values
        y_train = np.log10(train_df['mse_mean'].values)

        X_test  = test_df[feature_list].values
        y_test  = np.log10(test_df['mse_mean'].values)

        #do train/val split by chromosome 
        groups = train_df['CHROM'].values
        gkf = GroupKFold(n_splits=5)
        cv_folds = list(gkf.split(X_train, y_train, groups))

        enet_deg1_feat = Pipeline([
            ("poly", PolynomialFeatures(degree=1, interaction_only=True)),
            ("scaler", StandardScaler()),
            ("enet", ElasticNetCV(
                l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
                alphas=np.logspace(-4, 1, 30),
                cv=cv_folds, max_iter=50000,
                random_state=729, n_jobs=12, #verbose=1, 
            ))
        ])

        enet_deg1_feat.fit(X_train, y_train)

        best_model_deg1_feat = enet_deg1_feat.named_steps["enet"]
        best_alpha, best_l1=best_model_deg1_feat.alpha_, best_model_deg1_feat.l1_ratio_
        print("Best degree 1 alpha:", best_alpha)
        print("Best degree 1 l1_ratio:", best_l1)

        y_pred = enet_deg1_feat.predict(X_test)
        r2_test=r2_score(y_test, y_pred)
        mse_test=mean_squared_error(y_test, y_pred)

        print("degree 1 R2:", r2_test)
        print("degree 1 MSE:", mse_test)

        #add scores to df 
        final_df.append([feat_range, 1, best_alpha, best_l1, r2_test, mse_test])

        enet_cint_feat = Pipeline([
            ("poly", PolynomialFeatures(degree=2, interaction_only=True)),
            ("scaler", StandardScaler()),
            ("enet", ElasticNetCV(
                l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
                alphas=np.logspace(-4, 1, 30),
                cv=cv_folds, max_iter=50000,
                random_state=729, n_jobs=12, #verbose=1, 
            ))
        ])

        enet_cint_feat.fit(X_train, y_train)

        best_model_deg2_feat = enet_cint_feat.named_steps["enet"]
        best_alpha_2, best_l1_2=best_model_deg2_feat.alpha_, best_model_deg2_feat.l1_ratio_
        print("Best degree 2 alpha:", best_alpha_2)
        print("Best degree 2 l1_ratio:", best_l1_2)

        y_pred = enet_cint_feat.predict(X_test)
        r2_test_2=r2_score(y_test, y_pred)
        mse_test_2=mean_squared_error(y_test, y_pred)

        print("degree 2 R2:", r2_test_2)
        print("degree 2 MSE:", mse_test_2)

        #add scores to df 
        final_df.append([feat_range, 2, best_alpha_2, best_l1_2, r2_test_2, mse_test_2])

    ###############
    #combine final results 
    window_en_df=pd.DataFrame(final_df, columns=['window', 'interaction_degree', 'alpha', 'l1', 'r2', 'mse'])
    window_en_df.to_csv(f'{EN_dir}20260507_CVchromsplit_window_EN_performance_logMSE_{cell_type}.txt', sep='\t')

    #get feature names for 1Mb 
    poly1 = enet_deg1_feat.named_steps["poly"]
    enet1 = enet_deg1_feat.named_steps["enet"]
    feature_names = poly1.get_feature_names_out(feature_list)
    coef_df_deg1 = pd.DataFrame({"feature": feature_names,"coef": enet1.coef_}).sort_values("coef", key=np.abs, ascending=False)
    coef_df_deg1['degree']=1

    poly = enet_cint_feat.named_steps["poly"]
    enet = enet_cint_feat.named_steps["enet"]
    feature_names = poly.get_feature_names_out(feature_list)
    coef_df_deg2 = pd.DataFrame({"feature": feature_names,"coef": enet.coef_}).sort_values("coef", key=np.abs, ascending=False)
    coef_df_deg2['degree']=2

    coef_df_all=pd.concat([coef_df_deg1, coef_df_deg2], axis=0)
    coef_df_all.to_csv(f'{EN_dir}20260507_CVchromsplit_window_EN_coeffs_logMSE_{cell_type}.txt', sep='\t')

# ## Logistic Regression

#repeat for logistic reg
df_regress = no_na_alu

q99=np.quantile(df_regress['mse_mean'], 0.99)
q50=np.quantile(df_regress['mse_mean'], 0.50)

#make binary
alu_high=df_regress[df_regress['mse_mean']>=q99]
alu_low=df_regress[df_regress['mse_mean']<=q50]

alu_high['mse_binary']=1
alu_low['mse_binary']=0

alu_low_samp=alu_low.sample(n=len(alu_high), random_state=729)
alu_binary=pd.concat([alu_high, alu_low_samp], axis=0)

substrings_common = ['mse_binary', 'CHROM', 'orig_idx', 'GC_', 'mapp36', 'mapp100',
                     'phastCon100', 'alu_count', 'promoter_count',
                     'gene_PC_overlap', 'pseudogene_overlap', 'lincRNA_overlap',
                     'ZNF_overlap', 'repeats_overlap',
                     'phyloP100_highfrac', 'phyloP100_lowfrac',
                     'seg_dup_overlap', 'self_chain_overlap']
substrings_HFF = ['HFF_CTCF', 'HFF_H3K9me3', 'HFF_H3K27me3', 'IMR90_']
substrings_H1  = ['H1_CTCF',  'H1_H3K9me3',  'H1_H3K27me3',  'H1ESC_WGBS']

common_mask = alu_binary.columns.str.contains('|'.join(substrings_common), case=False)
hff_mask    = alu_binary.columns.str.contains('|'.join(substrings_HFF),    case=False)
h1_mask     = alu_binary.columns.str.contains('|'.join(substrings_H1),     case=False)

alu_binary_HFF = alu_binary.loc[:, common_mask | hff_mask]
alu_binary_H1  = alu_binary.loc[:, common_mask | h1_mask]

CT_dict_binary = {'HFF': alu_binary_HFF, 'H1': alu_binary_H1}

for cell_type in CT_dict_binary:
    df_regress = CT_dict_binary[cell_type]
    print(f'\n========== {cell_type} ==========')

    final_df_logreg = []
    roc_rows  = []   # long-format: one row per point on ROC curve
    pr_rows   = []   # long-format: one row per point on PR curve
    pred_rows = []   # raw (y_test, y_prob) pairs — most flexible for R
    coef_rows = []   # NEW: collect coefficients across all feat_range × degree

    for feat_range in ['alu', '1kb', '10kb', '100kb', '1Mb']:
        print(f'running {feat_range}...')
        cols_feat_range = df_regress.columns[df_regress.columns.str.endswith(feat_range)]
        feature_list = cols_feat_range.tolist()

        df_regress_nona = df_regress.copy()
        df_regress_nona[feature_list] = df_regress_nona[feature_list].fillna(0)

        test_chr = 'chr2'
        train_df = df_regress_nona[df_regress_nona['CHROM'] != test_chr].copy()
        test_df  = df_regress_nona[df_regress_nona['CHROM'] == test_chr].copy()

        X_train = train_df[feature_list].values
        y_train = train_df['mse_binary'].values
        X_test  = test_df[feature_list].values
        y_test  = test_df['mse_binary'].values

        print(f"  Train class balance: {y_train.mean():.2f} positive")
        print(f"  Test  class balance: {y_test.mean():.2f} positive")

        groups = train_df['CHROM'].values
        gkf = GroupKFold(n_splits=5)
        cv_folds = list(gkf.split(X_train, y_train, groups))

        for degree in [1, 2]:
            pipe = Pipeline([
                ("poly",   PolynomialFeatures(degree=degree, interaction_only=True)),
                ("scaler", StandardScaler()),
                ("clf",    LogisticRegressionCV(
                    Cs=np.logspace(-4, 1, 10),
                    cv=cv_folds,
                    penalty='elasticnet',
                    solver='saga',
                    l1_ratios=[0.1, 0.5, 0.9],
                    max_iter=20000,
                    random_state=729,
                    n_jobs=12,
                    scoring='roc_auc',
                    class_weight='balanced',
                ))
            ])
            pipe.fit(X_train, y_train)

            best_clf = pipe.named_steps["clf"]
            best_C   = best_clf.C_[0]
            best_l1  = best_clf.l1_ratio_[0]
            y_prob   = pipe.predict_proba(X_test)[:, 1]
            y_pred   = (y_prob >= 0.5).astype(int)
            auroc    = roc_auc_score(y_test, y_prob)
            auprc    = average_precision_score(y_test, y_prob)

            print(f"  Degree {degree} — C={best_C:.4f}, l1={best_l1:.2f}, "
                  f"AUROC={auroc:.3f}, AUPRC={auprc:.3f}")

            final_df_logreg.append([feat_range, degree, best_C, best_l1, auroc, auprc])

            # ROC curve points
            fpr, tpr, roc_thr = roc_curve(y_test, y_prob)
            for f, t, th in zip(fpr, tpr, roc_thr):
                roc_rows.append({
                    'cell_type': cell_type, 'feat_range': feat_range, 'degree': degree,
                    'fpr': f, 'tpr': t, 'threshold': th,
                })

            # PR curve points (precision_recall_curve returns one extra prec/rec
            # without a corresponding threshold — pad with NaN)
            prec, rec, pr_thr = precision_recall_curve(y_test, y_prob)
            pr_thr_padded = np.append(pr_thr, np.nan)
            for p, r, th in zip(prec, rec, pr_thr_padded):
                pr_rows.append({
                    'cell_type': cell_type, 'feat_range': feat_range, 'degree': degree,
                    'precision': p, 'recall': r, 'threshold': th,
                })

            # Raw predictions — most flexible if you want to recompute / bootstrap in R
            for yt, yp in zip(y_test, y_prob):
                pred_rows.append({
                    'cell_type': cell_type, 'feat_range': feat_range, 'degree': degree,
                    'y_test': int(yt), 'y_prob': float(yp),
                })

            #save coefficients 
            poly = pipe.named_steps["poly"]
            clf  = pipe.named_steps["clf"]
            feature_names = poly.get_feature_names_out(feature_list)
            for fname, c in zip(feature_names, clf.coef_[0]):
                coef_rows.append({
                    'cell_type':  cell_type,
                    'feat_range': feat_range,
                    'degree':     degree,
                    'feature':    fname,
                    'coef':       c,
                })

    #save summary stats 
    summary_df = pd.DataFrame(final_df_logreg, columns=['feat_range', 'degree', 'best_C', 'best_l1','auroc', 'auprc'])
    roc_df  = pd.DataFrame(roc_rows)
    pr_df   = pd.DataFrame(pr_rows)
    pred_df = pd.DataFrame(pred_rows)

    summary_df.to_csv(f'{EN_dir}20260507_logreg_summary_{cell_type}.tsv',sep='\t', index=False)
    roc_df.to_csv(f'{EN_dir}20260507_logreg_roc_curves_{cell_type}.txt', index=False)
    pr_df.to_csv(f'{EN_dir}20260507_logreg_pr_curves_{cell_type}.txt',  index=False)
    pred_df.to_csv(f'{EN_dir}20260507_logreg_predictions_{cell_type}.txt', index=False)

    #save coefficients 
    coef_df = pd.DataFrame(coef_rows)
    coef_df['abs_coef'] = coef_df['coef'].abs()
    coef_df = coef_df.sort_values(['feat_range', 'degree', 'abs_coef'],
                                  ascending=[True, True, False])
    coef_df.to_csv(f'{EN_dir}20260507_logreg_coefs_{cell_type}.txt',sep='\t', index=False)

final_df_logreg = []
for feat_range in ['alu', '1kb', '10kb', '100kb', '1Mb']:
#for feat_range in ['1Mb']:
    print(f'running {feat_range}...')
    cols_feat_range = df_regress.columns[df_regress.columns.str.endswith(feat_range)]
    
    feature_list = cols_feat_range.tolist()
    #feature_list=feature_list_subset

    df_regress_nona = alu_binary.copy()
    df_regress_nona[feature_list] = df_regress_nona[feature_list].fillna(0)

    test_chr = 'chr2'
    train_df = df_regress_nona[df_regress_nona['CHROM'] != test_chr].copy()
    test_df  = df_regress_nona[df_regress_nona['CHROM'] == test_chr].copy()

    X_train = train_df[feature_list].values
    y_train = train_df['mse_binary'].values

    X_test  = test_df[feature_list].values
    y_test  = test_df['mse_binary'].values

    print(f"  Train class balance: {y_train.mean():.2f} positive")
    print(f"  Test  class balance: {y_test.mean():.2f} positive")

    #do train/val split by chromosome 
    groups = train_df['CHROM'].values
    gkf = GroupKFold(n_splits=5)
    cv_folds = list(gkf.split(X_train, y_train, groups))

    #for degree in [1]:
    for degree in [1, 2]:
        interaction_only = True  # keeps same spirit as your enet setup
        pipe = Pipeline([
            ("poly",   PolynomialFeatures(degree=degree, interaction_only=interaction_only)),
            ("scaler", StandardScaler()),
            ("clf",    LogisticRegressionCV(
                #Cs=np.logspace(-4, 1, 30),       # inverse of alpha — 1/C
                Cs=np.logspace(-4, 1, 10),       # inverse of alpha — 1/C
                cv=cv_folds,
                penalty='elasticnet',
                solver='saga',                    # required for elasticnet
                #l1_ratios=[0.1, 0.3, 0.5, 0.7, 0.9],
                l1_ratios=[0.1, 0.5, 0.9],
                max_iter=20000,
                random_state=729,
                n_jobs=12,
                scoring='roc_auc',               # optimize for AUC
                class_weight='balanced',         # handles imbalanced classes
            ))
        ])

        pipe.fit(X_train, y_train)

        if degree == 1:
            lr_deg1_feat = pipe
        elif degree == 2:
            lr_deg2_feat = pipe

        best_clf = pipe.named_steps["clf"]
        best_C      = best_clf.C_[0]
        best_l1     = best_clf.l1_ratio_[0]
        print(f"  Degree {degree} — best C: {best_C:.4f}, best l1_ratio: {best_l1:.2f}")

        y_pred      = pipe.predict(X_test)
        y_prob      = pipe.predict_proba(X_test)[:, 1]

        auroc = roc_auc_score(y_test, y_prob)
        auprc = average_precision_score(y_test, y_prob)
        print(f"  Degree {degree} AUROC: {auroc:.3f}")
        print(f"  Degree {degree} AUPRC: {auprc:.3f}")
        print(classification_report(y_test, y_pred))

        final_df_logreg.append([feat_range, degree, best_C, best_l1, auroc, auprc])

final_df_logreg = pd.DataFrame(final_df_logreg, columns=['feat_range', 'degree', 'best_C', 'best_l1', 'auroc', 'auprc'])

window_logreg_df=pd.DataFrame(final_df_logreg, columns=['window', 'interaction_degree', 'best_C', 'l1', 'auroc', 'auprc'])
window_logreg_df.to_csv(f'{PROJECT_DIR}/results/20260212_elasticnet/20260317_CVchromsplit_window_logreg_performance_logMSE.txt', sep='\t')
window_logreg_df

# Degree 1
poly1 = lr_deg1_feat.named_steps["poly"]
lr1 = lr_deg1_feat.named_steps["clf"]
feature_names = poly1.get_feature_names_out(feature_list)
coef_df_deg1 = pd.DataFrame({"feature": feature_names, "coef": lr1.coef_[0]}).sort_values("coef", key=np.abs, ascending=False)
coef_df_deg1['degree'] = 1

# Degree 2
poly = lr_deg2_feat.named_steps["poly"]
lr2 = lr_deg2_feat.named_steps["clf"]
feature_names = poly.get_feature_names_out(feature_list)
coef_df_deg2 = pd.DataFrame({"feature": feature_names, "coef": lr2.coef_[0]}).sort_values("coef", key=np.abs, ascending=False)
coef_df_deg2['degree'] = 2

coef_df_all_logreg = pd.concat([coef_df_deg1, coef_df_deg2], axis=0)

coef_df_all_logreg.to_csv(f'{PROJECT_DIR}/results/20260212_elasticnet/20260317_CVchromsplit_window_logreg_coeffs.txt', sep='\t')

coef_df_deg2[0:20]

final_df_rf = []

for feat_range in ['alu', '1kb', '10kb', '100kb', '1Mb']:
    print(f'running {feat_range}...')
    cols_feat_range = df_regress.columns[df_regress.columns.str.endswith(feat_range)]

    feature_list = cols_feat_range.tolist()
    df_regress_nona = alu_binary.copy()
    df_regress_nona[feature_list] = df_regress_nona[feature_list].fillna(0)

    test_chr = 'chr2'
    train_df = df_regress_nona[df_regress_nona['CHROM'] != test_chr].copy()
    test_df  = df_regress_nona[df_regress_nona['CHROM'] == test_chr].copy()

    X_train = train_df[feature_list].values
    y_train = train_df['mse_binary'].values
    X_test  = test_df[feature_list].values
    y_test  = test_df['mse_binary'].values

    print(f"  Train class balance: {y_train.mean():.2f} positive")
    print(f"  Test  class balance: {y_test.mean():.2f} positive")

    groups = train_df['CHROM'].values
    gkf = GroupKFold(n_splits=5)
    cv_folds = list(gkf.split(X_train, y_train, groups))

    # Hyperparameter grid to search over
    param_grid = {
        'n_estimators': [300],
        'max_depth': [None, 10, 20],
        'max_features': ['sqrt', 0.3],
        'min_samples_leaf': [1, 5],
    }

    best_auroc = -1
    best_params = None
    best_model = None

    for n_est in param_grid['n_estimators']:
        for max_depth in param_grid['max_depth']:
            for max_feat in param_grid['max_features']:
                for min_leaf in param_grid['min_samples_leaf']:

                    clf = RandomForestClassifier(
                        n_estimators=n_est,
                        max_depth=max_depth,
                        max_features=max_feat,
                        min_samples_leaf=min_leaf,
                        class_weight='balanced',
                        random_state=729,
                        n_jobs=12,
                    )

                    # Manual GroupKFold CV to mirror logreg approach
                    fold_aurocs = []
                    for train_idx, val_idx in cv_folds:
                        clf.fit(X_train[train_idx], y_train[train_idx])
                        val_prob = clf.predict_proba(X_train[val_idx])[:, 1]
                        fold_aurocs.append(roc_auc_score(y_train[val_idx], val_prob))

                    mean_auroc = np.mean(fold_aurocs)

                    if mean_auroc > best_auroc:
                        best_auroc = mean_auroc
                        best_params = dict(
                            n_estimators=n_est,
                            max_depth=max_depth,
                            max_features=max_feat,
                            min_samples_leaf=min_leaf,
                        )

    print(f"  Best CV AUROC: {best_auroc:.3f} | Params: {best_params}")

    # Refit on full training set with best params
    best_model = RandomForestClassifier(
        **best_params,
        class_weight='balanced',
        random_state=729,
        n_jobs=12,
    )
    best_model.fit(X_train, y_train)

    y_pred = best_model.predict(X_test)
    y_prob = best_model.predict_proba(X_test)[:, 1]

    auroc = roc_auc_score(y_test, y_prob)
    auprc = average_precision_score(y_test, y_prob)

    print(f"  Test AUROC: {auroc:.3f}")
    print(f"  Test AUPRC: {auprc:.3f}")
    print(classification_report(y_test, y_pred))

    final_df_rf.append([
        feat_range,
        best_params['max_depth'],
        best_params['max_features'],
        best_params['min_samples_leaf'],
        best_auroc,
        auroc,
        auprc,
    ])

final_df_rf = pd.DataFrame(final_df_rf, columns=[
    'feat_range', 'max_depth', 'max_features', 'min_samples_leaf',
    'cv_auroc', 'test_auroc', 'test_auprc'
])

importances = best_model.feature_importances_
feat_imp_df = pd.DataFrame({
    'feature': feature_list,
    'importance': importances
}).sort_values('importance', ascending=False)
feat_imp_df

result = permutation_importance(
    best_model, X_test, y_test,
    n_repeats=10, random_state=729, n_jobs=12,
    scoring='roc_auc'
)
perm_imp_df = pd.DataFrame({
    'feature': feature_list,
    'importance_mean': result.importances_mean,
    'importance_std':  result.importances_std,
}).sort_values('importance_mean', ascending=False)

perm_imp_df

explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(X_test)  # shape: [n_samples, n_features, 2] for binary

# Summary plot — shows direction and spread
shap.summary_plot(shap_values[:, :, 1], X_test, feature_names=feature_list)

# Bar plot — mean absolute SHAP (like global importance)
shap.summary_plot(shap_values[:, :, 1], X_test, feature_names=feature_list, plot_type='bar')

rf = RandomForestClassifier(
    n_estimators=300,
    class_weight='balanced',
    random_state=729,
    n_jobs=12
)

selector = RFECV(
    estimator=rf,
    step=1,                    # remove 1 feature per iteration
    cv=cv_folds,               # your GroupKFold splits
    scoring='roc_auc',
    min_features_to_select=3,
    n_jobs=12
)
selector.fit(X_train, y_train)

print(f"Optimal number of features: {selector.n_features_}")
selected_features = [feature_list[i] for i in range(len(feature_list)) if selector.support_[i]]
print(selected_features)
