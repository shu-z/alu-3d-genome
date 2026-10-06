

import os
import os, psutil, io, gzip, time 
import warnings

from pybedtools import BedTool
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import pybedtools
import pysam
import seaborn as sns

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

os.chdir(f'{SUPREMO_DIR}/')
alu_dir=f'{PROJECT_DIR}/'
res_dir=f'{alu_dir}results/'
TF_dir=f'{alu_dir}results/20250514_TFenrich_RERUN/'
out_dir=f'{TF_dir}mutate/'
fa_path=f'{SUPREMO_DIR}/data/hg38.fa'
#tf_tracks=f'{POLLARD_DATA}/tf_binding/tfbs_motifs/JASPAR/UCSC_tracks/2024/hg38/'
#these are 2018 ones, it seems that I used these because the tracks are smaller 
tf_track_dir=f'{POLLARD_DATA}/tf_binding/tfbs_motifs/JASPAR/UCSC_tracks/2018/hg38/tsv/'

tracks=os.listdir(tf_track_dir) #579 tracks total 
print('lentracks', len(tracks))

alu_topbottom_path='202520513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt'
aludf=pd.read_csv(f'{res_dir}{alu_topbottom_path}', sep='\t', index_col=0)
aludf['POS']=aludf['POS'].astype(int)
aludf['END']=aludf['END'].astype(int)
alu_colorder = ['CHROM', 'POS', 'END']
aludf = aludf[alu_colorder + [col for col in aludf.columns if col not in alu_colorder]]

variant_df = aludf.copy()
TF_track_df=pd.DataFrame()

#run intersections 
start_time = time.time()
counter=0
for track_id in tracks:
    print(counter)

    #clear files from previous tracks 
    pybedtools.cleanup()

    #this is not ideal bc cleanup is removing it everytime, so need to keep loading it in
    alu_bedtool_unsort=BedTool.from_dataframe(aludf[alu_colorder])
    alu_bedtool=BedTool.sort(alu_bedtool_unsort)

    try:
        #MAKE SURE nrows is removed when actually running
        track_df=pd.read_csv(f'{tf_track_dir}{track_id}', sep='\t', header=None, skiprows=1,
                        names=['chr', 'start', 'end', 'tf', 'rel_score', '-log10(pval)', 'strand'])

        tf_name=track_df.iloc[1,3]
        print(tf_name)

        track_df['start']=track_df['start'].astype(int)
        track_df['end']=track_df['end'].astype(int)
        track_bedtool=BedTool.from_dataframe(track_df) 

        #intersect all alus 
        alu_intersect=alu_bedtool.intersect(track_bedtool, c=True, F=0.5).to_dataframe(header=None)

        #add number of intersections as a new row to dataframe 
        variant_df[f'{tf_name}']=alu_intersect['name']

        #create another dataframe with information about the track itself 
        motif_len=track_df['end'][0]-track_df['start'][0]
        TF_track_df = TF_track_df.append({'TF': tf_name, 
                                          'motif_len': motif_len, 
                                          'motif_occurence': len(track_df),  
                                          }, ignore_index=True)

    except:
        print(f'error in index {counter}')

#count up time for all intersections 
end_time = time.time()
elapsed_time = (end_time - start_time) / 60
print(f"Elapsed time: {elapsed_time} minutes")

variant_df.to_csv(f'{TF_dir}Alu_topbottom_TFoverlap_annot.txt', sep='\t')
TF_track_df.to_csv(f'{TF_dir}TF_tracks_info.txt', sep='\t')