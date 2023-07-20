configfile: "QC_config.yaml"

name=config["input_file"].keys()
output_path=config["output_path"]
ecc_master_path=config["ecc_master_path"]
threads=config["threads"]
paired=config["paired"]

import warnings

warnings.filterwarnings('ignore')

print(output_path)

def get_bwa_map_input_fastqs(wildcards):
    print(wildcards)
    return config["input_file"][wildcards.sample]

rule all:
    input:
        ##保留文件 这里要用expand 遍历
        expand('{output_path}/{name}',
               output_path=output_path,
               name=name)
##01.fastqc_trim
rule fastqc_trim:
    input:
        get_bwa_map_input_fastqs
        
    output: ##避免多个sample写入同一个output，必须添加路径
        output_path+'/{sample}/'
    
    log:
        output_path+'/logs/{sample}.01.log',
            
    threads: config["threads"]

    run:
        if paired:
            shell("./resource/QC/FastQC/fastqc {input} -t {threads} ; trim_galore \
            -q 30 -phred33 --stringency 3 --length 20 -e 0.1 --fastqc --paired {input} --gzip -o {output_path}/trim/")
        else:
            shell("./resource/QC/FastQC/fastqc {input} -t {threads} ; trim_galore \
            -q 30 -phred33 --stringency 3 --length 20 -e 0.1 --fastqc {input} --gzip -o {output_path}/trim/")


##02.picard_fastqc
#rule picard_fastqc:
    #input:
       # name='{sample}'_val_1.fq.gz'

    #params:

   # output:
        #output_path+''
    #run:
       ## if   '{input}_val_1.fq.gz' '{input}_val_2.fq.gz' else '{input}_val.fq.gz'

        










