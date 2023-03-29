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

def make_circlize_file(bed_file, geno, ecc_id):
    """
    bed_file = './data/result/cresil_new/cresil_result.analysis.bed'
    geno: 'hg38' or 'mm10'
    ecc_id: 21 int
    """
    path_share = '/'.join(bed_file.split('/')[:-1])
    for i in ['/05.circlize/']:
        _dir = path_share+i
        if not exists(_dir):
            makedirs(_dir)
            
    df = pd.read_csv(bed_file, sep='\t'
            , header=None)
    df.columns = ['Chr', 'Start', 'End', 'Count', 'eccID']
    df = df[df['eccID'].isin([ecc_id])]
    df.to_csv(path_share+'/05.circlize/circlize_ecc_peak.bed', sep='\t', header=None, index=None)
    
    ecc_bed_file = path_share+'/05.circlize/circlize_ecc_peak.bed'
    ref_gene_path = './resource/Analysis/reference/genes.10X.'+geno+'.bed'
    outfile_path= path_share+'/05.circlize/circlize_ecc_peak_gene.bed'
    
    os.system('bedtools intersect -a ' + ecc_bed_file + ' -b ' + ref_gene_path + ' -wa -wb -F 1 > ' + outfile_path)
    ## 计算 exon 完整在 ecc peak上的
    
    pdf_path = path_share+'/05.circlize/circlize_plot.pdf'
    os.system('Rscript ./analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path))
    print('Rscript ./analysis_code/circlize.R {0} {1} {2}'.format(ecc_bed_file, outfile_path, pdf_path))
    
    return None


