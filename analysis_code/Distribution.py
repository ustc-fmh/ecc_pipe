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
#mpl.rcParams['pdf.fonttype'] = 24
#mpl.rcParams["font.sans-serif"] = "Arial"
#plt.rcParams['font.size'] = 24
warnings.filterwarnings('ignore')

def plot_pie(df, font_size=9, output_path=None):
    """
    df: index: chr;  ['number']: number
    """
    ## rank
    df.columns = ['number']
    if len(df['number']) > 1:
#         df.index = df.index.astype('category')
        df_rank_list = df.index.values.sort_values(ascending=True, inplace=False)
        df = df.reindex(df_rank_list)
    else:
        print('only 1 chr be detect')
    ## plot
    fig = plt.figure(figsize = (8,8))
    ax = plt.subplot(1,1,1)
    df['number'].plot.pie(fontsize=font_size)
    ax.set_title('The proportion of chromosome types %',fontsize = font_size)
    
    if output_path != None:
        fig.savefig(output_path)

def premake_distribution(df):
    """
    df: result.columns in ['Length'], ['Count']
    """
    ## df_len_type_df
    bins = [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10]
    labels = ['$10^{' + str(len(str(bins[i]))-3) + '} - 10^{' + str(len(str(bins[i]))-2) + '}$' for i in range(len(bins)-1)]
    
    df.loc[:,'len_bin'] = pd.cut(df.Length, bins, labels = labels)
#     df_len_type_df = pd.DataFrame(df['Length'].value_counts()/ len(df))
    df_len_type_df = pd.DataFrame(df['Length'])
    df_len_type_df.columns = ['Ratio']
#     df_len_type_df.index.name='Length'
    df_len_type_df=df_len_type_df.reset_index()
    
    
    ## df_lenbin_type_df
    df_lenbin_type_df = pd.DataFrame(df['len_bin'].value_counts()/ len(df))
    df_lenbin_type_df.columns = ['Ratio']
    df_lenbin_type_df.index.name='Length'
    df_lenbin_type_df=df_lenbin_type_df.reset_index()
    
    ## df_len_count_df
#     df_group_len = df.groupby('Length').sum()
#     df_len_count_df = pd.DataFrame(df_group_len['Count'])
    df_len_count_df = pd.DataFrame(np.repeat(df['Length'],df['Count']))
#     df_len_count_df = pd.DataFrame(df_group_len['Count']/ sum(df_group_len.Count))
    df_len_count_df.columns = ['Ratio']
    df_len_count_df.index.name='Length'
    df_len_count_df=df_len_count_df.reset_index()
    
    ## df_lenbin_count_df
    df_group_lenbin = df.groupby('len_bin').sum()
    df_lenbin_count_df = pd.DataFrame(df_group_lenbin['Count']/sum(df_group_lenbin.Count))
    df_lenbin_count_df.columns = ['Ratio']
    df_lenbin_count_df.index.name='Length'
    df_lenbin_count_df=df_lenbin_count_df.reset_index()
    
    return df_len_type_df, df_lenbin_type_df, df_len_count_df, df_lenbin_count_df

def plot_distribution(len_df, lenbin_df, len_count_df, lenbin_count_df, name='default', output_path=None, xlim=2000):
    """
    len_df, lenbin_df
    name: str default 'cresil',
    output_path: None
    xlim: 2000
    _type: only plot distrubution for circlemap
    """
    
    fig = plt.figure(figsize=(25,20))
    fig = gs.GridSpec(10, 10, hspace=0.1, wspace=1)
    
    ax0 = plt.subplot(fig[0:3, 0:5])    
    ax1 = plt.subplot(fig[4:7, 0:5])
    ax2 = plt.subplot(fig[7:10, 0:5])
    ax3 = plt.subplot(fig[0:3, 5:10])
    ax4 = plt.subplot(fig[4:7, 5:10])
    ax5 = plt.subplot(fig[7:10, 5:10])
    
    ax0.set_title('Length Type Distribution', size=12, weight='semibold')
    ax1.set_title('Length Type Distribution (All)', size=12, weight='semibold')
    ax3.set_title('Length Count Distribution', size=12, weight='semibold')
    ax4.set_title('Length Count Distribution (All)', size=12, weight='semibold')
    
    
    sns.histplot(len_df, x='Ratio',
                 ax=ax0,
             log_scale=True,
             bins=[1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7],
             element="poly",
             thresh=10**7,alpha=0.6)
    sns.barplot(data=lenbin_df, x="Length", y='Ratio', errwidth=0.8,capsize=.2, ax=ax1, palette='Blues')
    sns.barplot(data=lenbin_df, x="Length", y='Ratio', errwidth=0.8,capsize=.2, ax=ax2, palette='Blues')
    
    sns.histplot(len_count_df, x='Ratio',
                 ax=ax3,
             log_scale=True,
             bins=[1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7],
             element="poly",
             thresh=10**7,alpha=0.6)
    sns.barplot(data=lenbin_count_df, x="Length", y='Ratio', errwidth=0.8,capsize=.2, ax=ax4, palette='Blues')
    sns.barplot(data=lenbin_count_df, x="Length", y='Ratio', errwidth=0.8,capsize=.2, ax=ax5, palette='Blues')
    
#     ax0.set_xlim([0, xlim])
    ax1.set_ylim([0.1, 0.7])
    ax2.set_ylim([0, 0.1])
#     ax3.set_xlim([0, xlim])
    ax4.set_ylim([0.1, 0.7])
    ax5.set_ylim([0, 0.1])
    
    ax1.set_xlabel('')
    ax1.set_xticks([])
    ax4.set_xlabel('')
    ax4.set_xticks([])
    
    ax1.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],labels=['0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','1.00'])
    ax4.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],labels=['0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','1.00'])
    
    if output_path != None:
        plt.savefig(output_path)
    
##后续homer一开始需要添加环境变量先
def annotate_homer(bed_path, geno):
    """
    bed_path:chr,start,end sep='\t'
    geno:mm10/hg38
    """
    
    output_path=bed_path.split('.bed')[0]+'.anno.tsv'
    _shell2call = 'annotatePeaks.pl {0} {1} > {2}'.format(bed_path, geno, output_path                                                                                                                                                 )
    sp.check_call(_shell2call, shell=True)
    
    return output_path

def get_anno_df(bed_path, geno):
    """
    bed_path:chr,start,end sep='\t'
    geno:mm10/hg38
    """
    anno_path = annotate_homer(bed_path, geno)
    anno_df = pd.read_csv(anno_path, sep='\t', index_col=0)
    ## drop na
    anno_df = anno_df.loc[anno_df["Annotation"].dropna().index, :]
    anno_df["Annotation_simple"] = anno_df["Detailed Annotation"]

    def _func(x):
        if x['Annotation_simple'] == x['Annotation']:
            if '(' in x['Annotation']:
                x['Annotation_simple'] = x['Annotation'].split('(')[0]
            else:
                x['Annotation_simple'] = x['Annotation']
        
#         if x['Annotation_simple'] == '':
#             x['Annotation_simple'] = x['Annotation'].split('(')[0]

        elif '|' in x['Annotation_simple']:
            x['Annotation_simple'] = x['Annotation_simple'].split('|')[1] ##取中间
            if '?' in x['Annotation_simple']:
                x['Annotation_simple'] = x['Annotation_simple'].split('?')[0]
                
        elif 'CpG-' in x['Annotation_simple']:
            x['Annotation_simple'] = 'CpG'
        else:
#             print(x['Annotation_simple'])
#             print(x['Annotation'])
            return x
        return x
    
    anno_df = anno_df.apply(func=_func, axis=1) 
    return pd.DataFrame(anno_df['Annotation_simple'].value_counts())

class distribution(object):
    def __init__(self,
                file_path: str,
                 _type:str,
                 geno:str,
                 bed_path=None,
                ):
        """
        file_path: AA,cresil,circlemap file path
        _type: str in ['AA', 'cresil', 'circlemap', 'other']
        geno: str in ['hg38', 'mm10']
        """
        self.save_path = '/'.join(file_path.split('/')[:-1])+'/ecc_pipe_result' ## the one up file path
        self.file_path = file_path
        self._type = _type
        self.geno = geno
        
        self.bed_path = bed_path
        
        self.print_deep = -1
        ## remake
        for i in ['/01.chr_distrbution/','/02.length_distrbution/',
                 '/03.homer_anno_distrbution/', '/04.db.annotation/']:
            _dir = self.save_path+i
            if not exists(_dir):
                makedirs(_dir)
        
        
    def deep_count(func):
        def func_wrapper(self, *args, **kwargs):
            self.print_deep += 1
            result = func(self, *args, **kwargs)
            self.print_deep -= 1
            return result
        return func_wrapper

    def myPrint(self, str):
        print('\t' * self.print_deep+str)
        
    @deep_count
    def make_bed_QC(self, _return=False):
        """
        _return: default:False, if True return file, bed_file
        """
        self.myPrint('Make chr bed & QC start!')
        
        if self._type == 'AA':
            self.file = pd.read_csv(self.file_path, sep='\t', index_col=0)
            self.file = self.file[self.file['Classification'].isin(['ecDNA'])]
            self.file = self.file[['Location', 'Captured interval length',
                                   'Feature median copy number', 'All genes', 'Feature BED file', 'AA amplicon number']]
            ## galaxy test update    
            ##self.file['Feature BED file'] = ['./data/hela_amplicon2_ecDNA_1_intervals.bed', './data/hela_amplicon4_ecDNA_1_intervals.bed']
            self.bed_file = pd.DataFrame()
            for path in self.file['Feature BED file']:
                middle = pd.read_csv(path, sep='\t', index_col=None, header=None)
                middle['Count'] = self.file[self.file['Feature BED file'].isin([path])]['Feature median copy number'].values
                middle['eccID'] = self.file[self.file['Feature BED file'].isin([path])]['AA amplicon number'].values
                ## AA eccid = amplicon
                self.bed_file = pd.concat([self.bed_file, middle], axis=0)
            self.bed_file.columns = ['Chr', 'Start', 'End', 'Count', 'eccID']
            
        elif self._type == 'circlemap':
            self.file = pd.read_csv(self.file_path,  sep='\t', header=None,
                              names=['Chr','Start','End','Discordants','Splits',
                                     'Score','Coverage','Std','Cov_Start','Cov_End','Coverage_Continuity'])
            self.file['Length'] = self.file['End'] - self.file['Start'] + 1
            self.file['Count'] = self.file['Discordants'] + self.file['Splits']

            ## QC https://github.com/iprada/Circle-Map/wiki/Circle-Map-Realign-output-files
            self.file = self.file[(self.file['Score']>50) & (self.file["Cov_Start"]>0.33) &
                 (self.file["Cov_End"]>0.33) &
                 (self.file["Coverage_Continuity"]<0.9) &
#                  (self.file["Coverage"]>0.5) & 
                 (self.file["Splits"]>2) & (self.file["Discordants"]>5) & (self.file['Length'] < 10000000)
                ]

#             self.file = self.file[(self.file['Score']>50) & 
#                                   (self.file["Cov_Start"]>0.33) &
#                  (self.file["Cov_End"]>0.33) &
#                  (self.file["Coverage_Continuity"]<0.5) &
#                  (self.file["Coverage"]>0.5) & 
#                  (self.file["Splits"]>2) & (self.file["Discordants"]>5)
#                 ]

            ##暂时无y,未考虑 mapping 其他一些，maybe bug
            self.file['Chr'] = self.file['Chr'].astype('category')
            self.file.Chr = self.file.Chr.cat.set_categories(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8'
                                                               ,'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16'
                                                               ,'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'])
            ## Chr Start End circlemap =
            self.file['eccID'] = self.file.index+1
            self.bed_file = self.file[['Chr', 'Start', 'End', 'Count', 'eccID']]
            
        elif self._type == 'cresil':
            self.file = pd.read_csv(self.file_path, sep='\t')
            self.file = self.file[self.file['eccdna_status'].isin(['cyclic'])]
            self.file['eccID'] = self.file.index+1
            
            cresil_chr_list = []
            cresil_length_list = []
            cresil_count_list = []
            cresil_id_list = []
            for i in self.file['merge_region']:
                count = self.file[self.file['merge_region'].isin([i])]['coverage'].values[0]
                eccID = self.file[self.file['merge_region'].isin([i])].index.values+1
                for j in i.split(','):
                    cresil_chr_list.append(j.split(':')[0])
                    cresil_length_list.append(j.split(':')[1])
                    cresil_count_list.append(count)
                    cresil_id_list.append(eccID[0])
            cresil_chr_df = pd.DataFrame([cresil_chr_list, cresil_length_list, cresil_count_list, cresil_id_list],
                                         index=['Chr', 'region', 'Count', 'eccID']).T
            cresil_chr_df['Chr'] = cresil_chr_df['Chr'].astype('category')
            cresil_chr_df.Chr = cresil_chr_df.Chr.cat.set_categories(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8'
                                                               ,'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16'
                                                               ,'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'])

            start_list = []
            end_list = []
            for i in cresil_chr_df['region']:
                #print(i.split('-'))
                for j in i.split('-'):
                    if '_' not in j and j!='':
                        start_list.append(j)
                    elif '_' in j:
                        end_list.append(j.split('_')[0])
                    else:
                        continue
            region_df = pd.DataFrame([start_list, end_list], index=['Start', 'End']).T
            cresil_chr_df = pd.concat([cresil_chr_df, region_df], axis=1)
            cresil_chr_df = cresil_chr_df.drop(['region'], axis=1)
            cresil_chr_df = cresil_chr_df[['Chr', 'Start', 'End', 'Count', 'eccID']]
            
            self.bed_file = cresil_chr_df
            
        elif self._type == 'other':
            if self.bed_path == None:
                print('Please set bed_path')
            else:
                self.file = pd.read_csv(self.file_path, sep='\t', index_col=0, header=0) 
                ## have header columns ['Chr', 'Start', 'End', 'Count', 'eccID', 'Length']
                self.bed_file = pd.read_csv(self.bed_path, sep='\t',header=None) 
                ## ['Chr', 'Start', 'End', 'Count', 'eccID'] not header columns
                self.bed_file.columns = ['Chr', 'Start', 'End', 'Count', 'eccID']
                self.bed_file['Chr'] = self.bed_file['Chr'].astype('category')
                self.bed_file.Chr = self.bed_file.Chr.cat.set_categories(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8'
                                                               ,'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16'
                                                               ,'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'])
        
        
        else:
            print("Please select str in ['AA', 'cresil', 'circlemap', 'other'] for _type paramater")
            
        
        ## save bed_file for homer
        self.bed_file.to_csv(self.save_path+'/'+self._type+'_result.analysis.bed',
                             sep='\t', header=False, index=False)
        self.file.to_csv(self.save_path+'/'+self._type+'_qc.txt', sep='\t', header=True,
                        index=True)
        
        ## for annotate_db
        self.bed_file_path =  self.save_path+'/'+self._type+'_result.analysis.bed'
        
        self.myPrint('Make chr bed & QC over!')
        
        if _return == 'False':
            return None
        else: 
            return self.file, self.bed_file
        
        
        
    @deep_count
    def plot_chr_distribution(self, _font_size=9, _return=False):
        """
        _font_size: default:9, font size for pie plot
        _return: default:False, if True return plot_chr_df
        """
        self.myPrint('Plot Chr distribution Start!')
        
        plot_chr_df = pd.DataFrame(self.bed_file['Chr'].value_counts())
        plot_pie(plot_chr_df, output_path=self.save_path+'/01.chr_distrbution/'+self._type+'_chr_distribution.pdf',
                font_size=_font_size)        
        plot_chr_df.to_csv(self.save_path+'/01.chr_distrbution/'+self._type+'_chr_distribution.csv',
                           header=True, index=True)
        
        self.myPrint('Plot Chr distribution over!')
        
        if _return == 'False':
            return None
        else: 
            return plot_chr_df
        
    @deep_count
    def plot_length_distribution(self, xlim=2000, _return=False):
        """
        xlim: distribution plot xlim(not bin) 
        _return: default:False, if True return len_df, lenbin_df, len_count_df, lenbin_count_df
        """
        self.myPrint('Plot Length distribution Start!')
        ## pre make
        if self._type=='circlemap':
            length_file = self.file
        elif self._type=='AA':
            length_file = self.file
            length_file['Length'] = length_file['Captured interval length']
            length_file['Count'] = length_file['Feature median copy number']
        elif self._type=='cresil':
            length_file = self.file
            length_file['Length'] = length_file['merge_len']
            length_file['Count'] = length_file['coverage']
        elif self._type=='other':
            length_file = self.file
        else:
            print("Please select str in ['AA', 'cresil', 'circlemap', 'other'] for _type paramater")
            
        len_df, lenbin_df, len_count_df, lenbin_count_df = premake_distribution(self.file)
        plot_distribution(len_df, lenbin_df,
                          len_count_df, lenbin_count_df,
                          name=self._type, xlim=xlim, 
                         output_path=self.save_path+'/02.length_distrbution/'+self._type+'_length_distribution.pdf')
        len_df.to_csv(self.save_path+'/02.length_distrbution/'+self._type+'_length_distribution.len.csv',
                     header=True, index=True)
        lenbin_df.to_csv(self.save_path+'/02.length_distrbution/'+self._type+'_length_distribution.lenbin.csv',
                     header=True, index=True)
        len_count_df.to_csv(self.save_path+'/02.length_distrbution/'+self._type+'_length_distribution.len_count.csv',
                     header=True, index=True)
        lenbin_count_df.to_csv(self.save_path+'/02.length_distrbution/'+self._type+'_length_distribution.lenbin_count.csv',
                     header=True, index=True)
        
        self.myPrint('Plot Length distribution over!')
        
        if _return == 'False':
            return None
        else: 
            return len_df, lenbin_df, len_count_df, lenbin_count_df
        
    @deep_count
    def plot_homer_anno_distribution(self, _return=False):
        """
        _return: default:False, if True return anno_df
        """
        self.myPrint('Plot homer anno distribution Start!')
        
        bed_path = self.save_path+'/'+self._type+'_result.analysis.bed'
        anno_df = get_anno_df(bed_path, self.geno)
        fig, ax = plt.subplots(figsize=(12,5))
        anno_df.plot(ax=ax, kind='bar', stacked=True, color=sns.color_palette("deep"))
        fig.savefig(self.save_path+'/03.homer_anno_distrbution/'+self._type+'_homer_anno_distribution.pdf')
        anno_df.to_csv(self.save_path+'/03.homer_anno_distrbution/'+self._type+'_homer_anno_distribution.csv',
                            header=True, index=True)
        
        self.myPrint('Plot homer anno distribution over!')
        
        if _return == 'False':
            return None
        else: 
            return anno_df
                
    @deep_count
    def annotate_db(self, trim=1):
        self.myPrint('Anno snp,se,e,eQTL Start!')
        
        if self.geno != 'hg38':
            print('This function just for hg38')
            
        ## annotate path
        snp_path = './resource/Analysis/reference/annotation/SNP.bed'
        se_path = './resource/Analysis/reference/annotation/superenhancer_dbsuper_sea_sedb.bed'
        e_path = './resource/Analysis/reference/annotation/enhancer_enhancerdb_sendb.bed'
        eQTL_path = './resource/Analysis/reference/annotation/eQTL.bed'
        
        snp_output_path = self.save_path+'/04.db.annotation/snp.result.bed'
        se_output_path = self.save_path+'/04.db.annotation/SuperEnhancer.result.bed'
        e_output_path = self.save_path+'/04.db.annotation/Enhancer.result.bed'
        eQTL_output_path = self.save_path+'/04.db.annotation/eQTL.result.bed'
        
        os.system('bedtools intersect -a ' + self.bed_file_path + ' -b ' + snp_path + ' -wa -wb -F '+ str(trim)+ ' > ' + snp_output_path)
        os.system('bedtools intersect -a ' + self.bed_file_path + ' -b ' + se_path + ' -wa -wb -F '+ str(trim)+ ' > ' + se_output_path)
        os.system('bedtools intersect -a ' + self.bed_file_path + ' -b ' + e_path + ' -wa -wb -F '+ str(trim)+ ' > ' + e_output_path)
        os.system('bedtools intersect -a ' + self.bed_file_path + ' -b ' + eQTL_path + ' -wa -wb -F '+ str(trim)+ ' > ' + eQTL_output_path)
        
        ## detect 0 
        if os.path.getsize(snp_output_path) == 0:
            snp_number = 0
        else:
            snp_df = pd.read_csv(snp_output_path, sep='\t', header=None)
            snp_number = snp_df.shape[0]
            
        if os.path.getsize(se_output_path) == 0:
            se_number = 0
        else:
            se_df = pd.read_csv(se_output_path, sep='\t', header=None)
            se_number = se_df.shape[0]
            
        if os.path.getsize(e_output_path) == 0:
            e_number = 0
        else:
            e_df = pd.read_csv(e_output_path, sep='\t', header=None)
            e_number = e_df.shape[0]
            
        if os.path.getsize(eQTL_output_path) == 0:
            eQTL_number = 0
        else:
            eQTL_df = pd.read_csv(eQTL_output_path, sep='\t', header=None)
            eQTL_number = eQTL_df.shape[0]
            
        ##plot number barplot    
        fig, ax = plt.subplots(figsize=(8,5))
        anno_df = pd.DataFrame([snp_number, se_number, e_number, eQTL_number],
                    index=['SNP', 'SuperEnhancer', 'Enhancer', 'eQTL'],
                    columns=['Annotation'])
        anno_df.plot(ax=ax, kind='bar', stacked=True, color=sns.color_palette("deep"))
        fig.savefig(self.save_path+'/04.db.annotation/Database_anno_number.pdf')
        anno_df.to_csv(self.save_path+'/04.db.annotation/Database_anno_number.csv',
                            header=True, index=True)
        
        
        self.myPrint('Anno snp,se,e,eQTL End!')
    
    @deep_count
    def run_fast(self, xlim=2000, trim=0.5):
        self.myPrint('Run fast Start!')
        self.make_bed_QC()
        self.plot_chr_distribution()
        self.plot_length_distribution(xlim)
        self.plot_homer_anno_distribution()
        if self.geno == 'hg38':
            self.annotate_db(trim=trim)
        
        self.myPrint('Run fast End!')
        
        

