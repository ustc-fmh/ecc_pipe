import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import matplotlib.gridspec as gs
import subprocess as sp
from os.path import join, exists
from os import makedirs
import seaborn as sns
import os
warnings.filterwarnings('ignore')

def make_circlize_file(bed_file, geno, ecc_id, ratio=0.5, _type='gene', ecc_pipe_path='None', R_env='None'):
    """
    bed_file = './data/result/cresil_new/cresil_result.analysis.bed'
    geno: 'hg38' or 'mm10'
    ecc_id: 21 int
    """
    path_share = '/'.join(bed_file.split('/')[:-1])
    for i in ['/06.circlize/']:
        _dir = path_share+i
        if not exists(_dir):
            makedirs(_dir)
            
    df = pd.read_csv(bed_file, sep='\t'
            , header=None)
    if len(df.columns) == 5:
        df.columns = ['Chr', 'Start', 'End', 'Count', 'eccID']
    elif len(df.columns) == 6:
        df.columns = ['Chr', 'Start', 'End', 'Count', 'eccID', 'Length']
    else:
        print("please set bed file have 5 cols in ['Chr', 'Start', 'End', 'Count', 'eccID']")
        
    df = df[df['eccID'].isin([ecc_id])]
    if df.shape[0] == 1:
        print('This eccDNA is simple, just 1 region')
    df.to_csv(path_share+'/06.circlize/circlize_ecc_peak.bed', sep='\t', header=None, index=None)
    
    ecc_bed_file = path_share+'/06.circlize/circlize_ecc_peak.bed'
    outfile_path= path_share+'/06.circlize/circlize_ecc_peak_gene.bed'
    
    ## set resource site
    if ecc_pipe_path=='None':
        ref_gene_path = './resource/Analysis/reference/genes.10X.'+geno+'.bed'
    else:
        ref_gene_path = ecc_pipe_path+'/resource/Analysis/reference/genes.10X.'+geno+'.bed'
    
    ### 
    if _type == 'gene':
        os.system('bedtools intersect -a ' + ecc_bed_file + ' -b ' + ref_gene_path + ' -wa -wb -F '+ str(ratio)+ ' > ' + outfile_path)
    elif _type == 'region':
        os.system('bedtools intersect -a ' + ecc_bed_file + ' -b ' + ref_gene_path + ' -wa -wb -f '+ str(ratio)+ ' > ' + outfile_path)
    else:
        print("Please set _type in ['gene', 'region']")
    
    ####计算 exon 完整在 ecc peak上的
    
    pdf_path = path_share+'/06.circlize/circlize_plot.pdf'
    
    if ecc_pipe_path=='None':
        if R_env == 'None':
            shell_code = 'Rscript ./analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path)
        else:
            shell_code = R_env+' ./analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path)
    else:
        if R_env == 'None':
            shell_code = 'Rscript '+ecc_pipe_path+'/analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path)
        else:
            shell_code = R_env+' '+ecc_pipe_path+'/analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path)
            
    os.system(shell_code)
    #print(shell_code)
    
    return None


