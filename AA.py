configfile: "AA_config.yaml"

##define params
if config['user_ref'] == 'None':
    if config['reference']=='hg38':
        reference = 'GRCh38'
    else:
        reference = config['reference']
else:
    print("AA don't accept other species; Please set user_ref == 'None'")
    
name=config["input_file"].keys()
output_path=config["output_path"]
ecc_master_path=config["ecc_master_path"]
threads=config["threads"]

cnvkit_dir=config["cnvkit_dir"]
rscript=config["rscript"]
python_path=config["python_path"]
cngain=config["cngain"]
cnsize=config["cnsize"]

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
##01.AA
rule AA:
    input:
        get_bwa_map_input_fastqs

    params:
        name='{sample}',
        ref=reference,
        cnvkit_dir=cnvkit_dir,
        rscript=rscript,
        python_path=python_path,
        cngain=cngain,
        cnsize=cnsize,
        
    output: ##避免多个sample写入同一个output，必须添加路径
        output_path+'/{sample}/'
    
    log:
        output_path+'/logs/{sample}.01.log',
            
    threads: config["threads"]
    shell:
        "python ./resource/AA/AmpliconSuite-pipeline/PrepareAA.py \
            -o {output} -s {params.name} -t {threads} \
            --cnvkit_dir {params.cnvkit_dir} \
            --fastqs {input} --ref {params.ref} \
            --cngain {params.cngain} --cnsize_min {params.cnsize} \
            --rscript_path {params.rscript} --python3_path {params.python_path} \
            --run_AA --run_AC"
        


        










