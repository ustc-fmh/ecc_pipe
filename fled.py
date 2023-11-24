configfile: "fled_config.yaml"
    
##define params
name=config["input_file"].keys()
output_path=config["output_path"]
ecc_master_path=config["ecc_master_path"]
threads=config["threads"]

##take cresil ref
if config['user_ref'] == 'None':
    config_ref = config['reference']
    reference=ecc_master_path+'/resource/cresil-master/reference/'+config_ref+'/'+config_ref+'.fa'
    
else:
    config_ref = config['user_ref']
    reference=ecc_master_path+'/resource/user_add/reference/'+config_ref+'.fa'

import warnings
warnings.filterwarnings('ignore')

def get_bwa_map_input_fastqs(wildcards):
    print(wildcards)
    return config["input_file"][wildcards.sample]

rule all:
    input:
        ##保留文件 这里要用expand 遍历
        expand("{output_path}/01.mapped_reads/{name}_sort.bam",
               output_path=output_path,
               name=name),
        expand("{output_path}/02.Circle_detection/{name}/",
               output_path=output_path,
               name=name),

##01. minimap2 mapping
rule minimap2_map:
    input:
        expand("{reference}", reference=reference),
        get_bwa_map_input_fastqs
        
    output:
        output_path+'/01.mapped_reads/{sample}.sam',
    
    log:
        output_path+'/logs/{sample}.01.log',
            
    threads: config["threads"]
    shell:
        "minimap2 -ax map-ont \
                -t {threads} \
            {input} > {output}"

##02.samtools_sort
rule samtools_sort:
    input:
        output_path+'/01.mapped_reads/{sample}.sam',
    output:
        output_path+'/01.mapped_reads/{sample}_sort.bam',
    log:
        output_path+'/logs/{sample}.02.log',

    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} -o {output} {input} ; "
        "samtools index -@ {threads} {output}"

##03.FLED_Detection
rule FLED_Detection:
    input:
        bam = output_path+'/01.mapped_reads/{sample}_sort.bam',
        fastq = get_bwa_map_input_fastqs
    output:
        path = output_path+'/02.Circle_detection/{sample}/'
    log:
        output_path+'/logs/{sample}.03.log',

    threads: config["threads"]
    shell:
        "FLED Detection -i {input.bam} -fq {input.fastq} -o circle_detected_ -dir {output.path} -t {threads} "










