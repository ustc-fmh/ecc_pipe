configfile: "cresil_config.yaml"
    
##define params   
name=config["input_file"].keys()
output_path=config["output_path"]
ecc_master_path=config["ecc_master_path"]
threads=config["threads"]
wgs=config["wgs"]

## user add ref
if config['user_ref'] == 'None':
    config_ref = config['reference']
    ref=ecc_master_path+'/resource/cresil-master/reference/'+config_ref+'/'+config_ref+'.mmi'
    fa=ecc_master_path+'/resource/cresil-master/reference/'+config_ref+'/'+config_ref+'.fa'
    fai=ecc_master_path+'/resource/cresil-master/reference/'+config_ref+'/'+config_ref+'.fa.fai'
    
else:
    config_ref = config['user_ref']
    ref=ecc_master_path+'/resource/user_add/reference/'+config_ref+'.mmi'
    fa=ecc_master_path+'/resource/user_add/reference/'+config_ref+'.fa'
    fai=ecc_master_path+'/resource/user_add/reference/'+config_ref+'.fa.fai'

#rp=ecc_master_path+'/resource/cresil-master/reference/'+config['reference']+'/reference.rmsk.bed'
#cg=ecc_master_path+'/resource/cresil-master/reference/'+config['reference']+'/reference.cpg.bed'
#gb=ecc_master_path+'/resource/cresil-master/reference/'+config['reference']+'/reference.gene.bed'
## cresil annotate -t {threads} -rp {params.rp} -cg {params.cg} -gb {params.gb} -identify {output}/eccDNA_final.txt
## rp=rp,cg=cg,gb=gb,

import warnings
warnings.filterwarnings('ignore')


def get_bwa_map_input_fastqs(wildcards):
    print(wildcards)
    return config["input_file"][wildcards.sample]

rule all:
    input:
        ##保留文件 这里要用expand 遍历
        expand('{output_path}/{name}',
               output_path=output_path,
               name=name)
##01.cresil
rule cresil:
    input:
        get_bwa_map_input_fastqs

    params:
        ref=ref,
        fa=fa,
        fai=fai,
        
    output: ##避免多个sample写入同一个output，必须添加路径
        output_path+'/{sample}/'
    
    log:
        output_path+'/logs/{sample}.01.log',
            
    threads: config["threads"]
    run:
        if wgs:
            shell("cresil trim -t {threads} -fq {input} -r {params.ref} -o {output};cresil identify_wgls -t {threads} \
            -r {params.ref} -fa {params.fa} -fai {params.fai} -fq {input} -trim {output}/trim.txt")
        else:
            shell("cresil trim -t {threads} -fq {input} -r {params.ref} -o {output};cresil identify -t {threads} \
            -fa {params.fa} -fai {params.fai} -fq {input} -trim {output}/trim.txt")
        
        


        










