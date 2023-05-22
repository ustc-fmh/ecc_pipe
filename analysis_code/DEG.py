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

import math
def volcano_plot(df_DEG, pvalue=0.05, log2FC=1, xlim=10, ylim=5, save=None, anno=True):
    
    result = pd.DataFrame()
    result['x'] = df_DEG['log2FoldChange']
    result['y'] = -df_DEG['pvalue'].apply(np.log10)

    #分组为up, normal, down
    result['group'] = 'black'
    result.loc[(result.x > log2FC) & (result.y > -math.log10(pvalue)), 'group'] = 'tab:red' #x=-+pvalue直接截断
    result.loc[(result.x < -log2FC) & (result.y > -math.log10(pvalue)), 'group'] = 'tab:blue' #x=-+pvalue直接截断
    result.loc[result.y < -math.log10(pvalue), 'group'] = 'dimgrey' #阈值以下点为灰色

    print(result[result['group'] == 'tab:red'].shape)
    print(result[result['group'] == 'tab:blue'].shape)

    #确定坐标轴显示范围
    xmin=-xlim
    xmax=xlim
    ymin=0
    ymax=ylim

    #绘制散点图
#     sns.set_style("whitegrid", {'axes.grid' : False})
    fig = plt.figure(figsize=plt.figaspect(5/6)) #确定fig比例（h/w）
    ax = fig.add_subplot()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='DEG')
    result_color = result[result['group'].isin(['tab:red', 'tab:blue'])]
    result_grey = result[result['group'].isin(['dimgrey'])]


    ax.scatter(result_grey['x'], result_grey['y'], s=10, c=result_grey['group'], alpha=0.05)
    ax.scatter(result_color['x'], result_color['y'], s=10, c=result_color['group'])
    
    ax.set_ylabel('-Log10(P-value)')
    ax.set_xlabel('Log2(FoldChange)')
    ax.spines['right'].set_visible(False) #去掉右边框
    ax.spines['top'].set_visible(False) #去掉上边框
    
    if anno != False:
        for xyz in zip(result_color['x'], result_color['y'], result_color.index.values):
            xy = (xyz[0], xyz[1])
            plt.annotate(xyz[2], xy=xy, xytext=(0, 0), textcoords='offset points', size=10)
    
    
    #水平和竖直线
    ax.vlines(-log2FC, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖直线
    ax.vlines(log2FC, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖直线
    ax.hlines(-math.log10(pvalue), xmin, xmax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖水平线

    plt.show()
    
    if save != None:
        fig.savefig(save, dpi=300, bbox_inches='tight')

def make_ecc_gene_martix(geno, bed_file_path, trim=1, _type='gene'):
    """
    geno: str in ['hg38', 'mm10']
    bed_file_path: peak list bed file path and Count
    trim: bedtools intersect ratio, 0-1
    _type: all gene in ecc region or ecc region in gene / str in ['gene', 'region']
    """
    ## 已修改
    ref_gene_path = './resource/Analysis/reference/genes.10X.'+geno+'.all.bed'

    outfile_up_path = '/'.join(bed_file_path.split('/')[:-1])
    outfile_path = outfile_up_path+'/ecc_gene_number.bed'
    
    if _type == 'gene':
        os.system('bedtools intersect -a ' + bed_file_path + ' -b ' + ref_gene_path + ' -wa -wb -F '+ str(trim)+ ' > ' + outfile_path)
    elif _type == 'region':
        os.system('bedtools intersect -a ' + bed_file_path + ' -b ' + ref_gene_path + ' -wa -wb -f '+ str(trim)+ ' > ' + outfile_path)
    else:
        print("Please set _type in ['gene', 'region']")
        
    bed_file = pd.read_csv(outfile_path, sep='\t', header=None)
    bed_file.columns=['Chr_peak', 'Start_peak', 'End_peak', 'Count', 'eccID', 'Chr_gene', 'Start_gene', 'End_gene', 'gene']
    
    
    gene_number_df = pd.DataFrame(bed_file.groupby(['gene'])['Count'].sum())
    gene_number_df.columns = ['gene_number']
    gene_number_df.to_csv(outfile_up_path+'/ecc_gene_number.csv', header=True, index=True)
    
    return gene_number_df

def plot_go_kegg_clusterprofile(result_path, top_number=10):
    """
    result_path: 'data/result/deg_test/02.ecc_deg_go_kegg'
    """
    for name in os.listdir(result_path):
        if name.split('.')[-1] == 'csv':
            #print(name)
            plot_path=result_path+'/'+name.split('.')[0]+'.pdf'
            
            df = pd.read_csv(result_path+'/'+name)
            if df.shape[0] < top_number:
                top_number=df.shape[0]
            df = df[:top_number]
            df = df[['Description', 'p.adjust']]
            df['-log10_p.adjust'] =  -df['p.adjust'].apply(np.log10)
            #print(df)
            #print(plot_path)
            ## init
            fig = plt.figure()
            sns.barplot(x= '-log10_p.adjust',y= "Description", data = df, palette="Reds_d")
            plt.savefig(plot_path, bbox_inches='tight')
            plt.show()
        
class ecc_gene_number_deg(object):
    def __init__(self,
                path_share: str,
                 group_file_path: str,
                 geno:str,
                 trim=1,
                 _type='gene',
                ):
        """
        path_share: './data/result/deg_test/'
        group_file_path: './data/result/deg_test/group.txt'
        geno: str in ['hg38', 'mm10']
        trim: bedtools intersect ratio, 0-1
        _type: all gene in ecc region or ecc region in gene / str in ['gene', 'region']
        """
        self.path_share = path_share ## the one up file path
        self.group_file_path = group_file_path
        self.group_file = pd.read_csv(group_file_path,
                                      sep='\t', header=None)
        self.group_file.columns=['name', 'group', 'tool']
        
        self.trim = trim
        self._type = _type
        self.geno = geno
        
        self.print_deep = -1
        ## remake
        for i in ['/01.ecc_deg_output/', '/02.ecc_deg_go_kegg/']:
            _dir = self.path_share+i
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
    def make_ecc_number_matrix(self):
        self.myPrint('Make ecc gene number matrix start!')
        
        for _name in self.group_file['group'].values:
            tool = self.group_file[self.group_file['group'].isin([_name])]['tool'].values[0]
            bed_file_path = self.path_share + '/'+_name + '/ecc_pipe_result/'+tool+'_result.analysis.bed'
            make_ecc_gene_martix(self.geno, bed_file_path,
                             trim=self.trim,
                            _type=self._type) ## make matrix
        
        ##merge
        df_merge = pd.DataFrame()
        for _name in self.group_file['group'].values:
            matrix_path = self.path_share + '/' +_name + '/ecc_pipe_result/ecc_gene_number.csv'
            df_middle = pd.read_csv(matrix_path, index_col=0, header=0)
            df_middle.columns = [_name]
            df_merge = pd.concat([df_merge, df_middle], axis=1)
        df_merge = df_merge.fillna(0)
        df_merge = df_merge.applymap(lambda x: int(x))
        
        df_merge.to_csv(self.path_share+'/01.ecc_deg_output/count.txt', sep='\t', header=True, index=True)
        
        self.count_file_path = self.path_share+'/01.ecc_deg_output/count.txt'
        self.myPrint('Make ecc gene number matrix end!')
        
    @deep_count
    def deseq2_run(self, pvalue=0.05, log2fc=1, xlim=10, ylim=5, anno=True):
        """
        xxx
        """
        self.myPrint('deseq2 run start!')
        outputdir_path = self.path_share+'/01.ecc_deg_output/'
        os.system('Rscript ./analysis_code/deseq2.R {0} {1} {2}'.format(self.count_file_path,
                                                          self.group_file_path,
                                                          outputdir_path))
        print('Rscript ./analysis_code/deseq2.R {0} {1} {2}'.format(self.count_file_path,
                                                          self.group_file_path,
                                                          outputdir_path))
        ## 可补充一个python 火山图绘制
        deg_df_path = outputdir_path+'/deseq2_result.csv'
        deg_df = pd.read_csv(deg_df_path, index_col=0)
        volcano_path = outputdir_path+'/deseq2.volcano.pdf'
        
        volcano_plot(deg_df, pvalue=pvalue, log2FC=log2fc, xlim=xlim, ylim=ylim,
             save=volcano_path,
            anno=anno)
        
        self.myPrint('deseq2 run end!')
        
    def clusterprofile_run(self, pvalue=0.05, log2fc=1):
        """
        xxx
        """
        self.myPrint('clusterprofile run start!')
        ###cluster profile
        deg_result = self.path_share+'/01.ecc_deg_output/deseq2_result.csv'
        outputdir_path = self.path_share+'/02.ecc_deg_go_kegg/'
        
        os.system('Rscript ./analysis_code/clusterprofile.R {0} {1} {2} {3} {4}'.format(deg_result,
                                                                             outputdir_path,
                                                                             pvalue,
                                                                             log2fc,
                                                                             self.geno))
        
        ## 增加 python 绘图 barplot
        print('Rscript ./analysis_code/clusterprofile.R {0} {1} {2} {3} {4}'.format(deg_result,
                                                                             outputdir_path,
                                                                             pvalue,
                                                                             log2fc,
                                                                             self.geno))
        
        plot_go_kegg_clusterprofile(outputdir_path)
        self.myPrint('clusterprofile run end!')
        
        
    @deep_count
    def run_fast(self):
        self.myPrint('Run fast Start!')
        
        self.make_ecc_number_matrix()
        self.deseq2_run()
        self.clusterprofile_run()

        self.myPrint('Run fast End!')
        

