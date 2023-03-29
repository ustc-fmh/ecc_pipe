import os
import sys
from optparse import OptionParser
from optparse import OptionGroup

sys.path.append('./analysis_code')
from Distribution import *
from DEG import *
from Circlize import *


#######  Options  #######
usage="usage: %prog [options][inputs]"
parser=OptionParser(usage=usage, version="%prog 1.0")

## general params
parser.add_option("-n",type='int',help="Set the max thread") ## for QC/Detect

## QC params
parser.add_option("--QC",action='store_true',default=False,
                      help="use fastqc,trim_glores,picards")

## Detect and Analysis share params
parser.add_option("--tool",type='string',help="set tools to analysis eccDNA")

## Detect params
parser.add_option("--Detect",action='store_true',default=False,
                      help="use tools for detect eccDNA")

## Analysis params
parser.add_option("--bed_file",type='string',help="*_result.analysis.bed")
parser.add_option("--ecc_id",type='int',help="ecc id (int) in the bed_file")
## Circlize

parser.add_option("--geno",type='string',help="hg38 or mm10")
parser.add_option("--trim",type='string',help="bedtools intersect ratio (0-1)")
parser.add_option("--path_share",type='string',help="DEG input file share path")
parser.add_option("--group_file",type='string',help="txt file whcih contain group info")
parser.add_option("--count_type",type='string',help="gene or region")
## DEG

parser.add_option("--peak_path",type='string',help="eccDNA result peak file path") ## tool = other
parser.add_option("--file_path",type='string',help="eccDNA result file path")
## Distribution

parser.add_option("--mode",type='string',help="str in ['Distrbution', 'DEG', 'Visualize'] ")
parser.add_option("--Analysis",action='store_true',default=False,
                      help="plot & analysis")

(options, args) = parser.parse_args()

#print(type(options))
#print(options)
#print(options.tool)
def Detect():
    """
    Detect
    """
    tool_set = options.tool
    threads_number = options.n

    if tool_set == 'circlemap':
        print('run circle-map for eccDNA')
        #os.system('snakemake --cores {0} -n -s circlemap.py')
        #os.system('snakemake --cores {0} -s circlemap.py')
        os.system('snakemake --cores {0} -n -s circlemap.py'.format(threads_number)) ## print log  test run
        os.system('snakemake --unlock -s circlemap.py') ## unlock
        os.system('snakemake --cores {0} -s circlemap.py'.format(threads_number)) ## --report --all-temp
    
    elif tool_set == 'AA':
        print('run AA for eccDNA')
        os.system('snakemake --cores {0} -n -s AA.py'.format(threads_number)) ## print log 试运行
        os.system('snakemake --unlock -s AA.py') ## unlock
        os.system('snakemake --cores {0} -s AA.py'.format(threads_number)) ## --report --all-temp

        ##detect index.html in pwd and rm
        if os.path.exists("./index.html"):
            os.remove('./index.html')

    
    elif tool_set == 'cresil':
        print('run cresil for eccDNA')
        os.system('snakemake --cores {0} -n -s cresil.py'.format(threads_number)) ## print log test run
        os.system('snakemake --unlock -s cresil.py') ## unlock
        os.system('snakemake --cores {0} -s cresil.py'.format(threads_number)) ## --report --all-temp
    
    else:
        print("Please set tool in ['circlemap', 'AA', 'cresil']")
        

def QC():
    """
    fastqc/trim_glores/picards
    """
    threads_number = options.n
    print('run QC')
    os.system('snakemake --cores {0} -n -s QC.py'.format(threads_number)) ## print log  test run
    os.system('snakemake --unlock -s QC.py') ## unlock
    os.system('snakemake --cores {0} -s QC.py'.format(threads_number)) ## --report --all-temp


def Analysis():
    """
    analysis & plot result
    """
    mode = options.mode
    if mode == 'Distribution':
        geno = options.geno
        file_path = options.file_path
        _type = options.tool
        print(geno)
        print(file_path)
        print(_type)
        if _type == 'other':
            peak_path = options.peak_path
            _distribution = distribution(file_path=file_path, _type=_type, geno=geno, bed_path=peak_path)
        elif _type in ['circlemap', 'AA', 'cresil']:
            _distribution = distribution(file_path=file_path, _type=_type, geno=geno)
        else:
            print("Please set tool in ['circlemap', 'AA', 'cresil', 'other']")
    
        trim = options.trim
        _distribution.run_fast(trim=trim)

        
    elif mode == 'DEG':
        #print('DEG')
        path_share = options.path_share
        group_file = options.group_file
        geno = options.geno
        trim = options.trim
        _type = options.count_type
        if _type in ['gene', 'region']:
            deg_class = ecc_gene_number_deg(path_share=path_share,
                                    group_file_path=group_file,
                                    geno=geno,
                                    trim=trim,
                                   _type=_type)
            deg_class.run_fast()
        else:
            print("Please set count_type in ['gene', 'region']")
        
        
    elif mode == 'Visualize':
        #print('Circlize')
        geno = options.geno
        bed_file = options.bed_file
        ecc_id = options.ecc_id
        make_circlize_file(bed_file, geno, ecc_id)
        
    else:
        print("Please set mode in ['Distribution', 'DEG', 'Visualize']")	
    

 

def main():
    if options.QC: QC()
    if options.Detect: Detect()
    if options.Analysis: Analysis()
    if(options.QC==False) & (options.Detect==False) & (options.Analysis==False):
        print("Please select a function in ['--QC', '--Detect', '--Analysis']")

main()

    
