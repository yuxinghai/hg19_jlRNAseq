import os
project="/hg19_jlRNAseq"
base="/data2/zhoulab/yuxinghai"
proj_dir= base+project
work_dir=proj_dir+"/results"
anno_dir= proj_dir+ "/anno"
bw_dir=work_dir+"/bw"
raw_dir="/zhounas/zhoulab/backup/xhyu/20171012JL/HB160055-01武汉大学12个人lncRNA建库测序过滤任务单_VIP/ANHB160055_HB160055-01_2017-01-09/Cleandata"
script_dir=proj_dir+"/script"
results_dir=work_dir+"/03_featurecount"
samples = os.listdir(raw_dir)
samples.sort()
py27_dir = "/data1/zhoulab/yuxinghai/.conda/envs/py27/bin"
gene_list=["up","down","sig"]
period = []
chrom_size= anno_dir + "/hg19.chrom.sizes"


rule all:
    input:
        expand("{anno_dir}/rRNA/SAindex",anno_dir=anno_dir),
        expand("{anno_dir}/hg19_index/SAindex",anno_dir=anno_dir),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort", work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam",work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam",work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam",work_dir=work_dir, dataset=samples),
        expand("{work_dir}/03_featurecount/F1_featureCounts.txt", work_dir=work_dir),
        expand("{work_dir}/03_featurecount/F1_featureStat.log", work_dir=work_dir),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam.bai",work_dir=work_dir,dataset=samples),
        expand("{work_dir}/bw/{dataset}.Forward.bw",work_dir=work_dir,dataset=samples),
        expand("{results_dir}/figure/distance_PNG_300_lengend_dpi.png",results_dir=results_dir),
        expand("{results_dir}/DEgene/PSF/ano_downgene.tsv",results_dir=results_dir),
        expand("{results_dir}/DEgene/NONO/ano_downgene.tsv",results_dir=results_dir),
        expand("{results_dir}/DEgene/PSPC1/ano_downgene.tsv",results_dir=results_dir),
        expand("{results_dir}/DEgene/NEAT1_V2/ano_downgene.tsv",results_dir=results_dir),
        expand("{results_dir}/DEgene/NEAT1_NT/ano_downgene.tsv",results_dir=results_dir),





rule build_index:
    input:
        rRNA_fa="{anno_dir}/hg19_encodev19_rRNA.fa",
        gtf="{anno_dir}/gencode.v19.annotation.gtf",
        fa="{anno_dir}/hg19.fa"

    output:
        rRNA_index="{anno_dir}/rRNA/SAindex",
        genome_index="{anno_dir}/hg19_index/SAindex"

    threads: 20

    shell: """
        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {anno_dir}/rRNA --genomeFastaFiles {input.rRNA_fa} --genomeSAindexNbases 6 --outFileNamePrefix {anno_dir}/rRNA/hg19_rRNA

        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {anno_dir}/hg19_index --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --genomeSAindexNbases 12 --outFileNamePrefix {anno_dir}/hg19_index/hg19
        """

rule rRNA_map:
    input:
        R1 = raw_dir+"/{dataset}/{dataset}_R1.fq.gz",
        R2 = raw_dir+"/{dataset}/{dataset}_R2.fq.gz"

    output:
        "{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1",
        "{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2"

    threads: 20

    shell:
        """STAR --genomeDir {anno_dir}/rRNA --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c --runThreadN 20 --outFileNamePrefix {work_dir}/02_mapping/rRNA/{wildcards.dataset}/{wildcards.dataset} --outReadsUnmapped Fastx --outFilterMatchNmin 40 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
        """


rule sort_STAR_result:
    input:
        un1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1",
        un2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2"


    output:
        fastq1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort",
        fastq2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort"

    shell:
        """
        cat {input.un1}| paste -d " " - - - - | sort -k1,1 -S 10G | tr ' ' '\n' > {output.fastq1}
        cat {input.un2}| paste -d " " - - - - | sort -k1,1 -S 10G | tr ' ' '\n' > {output.fastq2}
        """

rule genome_map:
    input:
        f_fastq1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort",
        f_fastq2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort"
    output:
        f_Sam="{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam"

    shell:
        """
        STAR --genomeDir {anno_dir}/hg19_index --readFilesIn {input.f_fastq1} {input.f_fastq2} --runThreadN 20 --outFileNamePrefix {work_dir}/02_mapping/sorted/{wildcards.dataset}/{wildcards.dataset}_Sort --outReadsUnmapped Fastx --outFilterMatchNmin 40 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
        """

rule sam_2_uniquebam:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam"

    output:
        u_bam="{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam"

    shell:
        "samtools view -@ 8 -h -bq 255 {input} > {output.u_bam}"

rule uniquebam_sort_and_index:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam"

    output:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam"

    shell:
        "samtools sort -@ 20 {input} -o {output}"


rule bam_index_and_2bw:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam"

    output:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam.bai",
        "{work_dir}/bw/{dataset}.Forward.bw",
        "{work_dir}/bw/{dataset}.Reverse.bw"

    shell:
        """
        samtools index {input}
        {py27_dir}/bam2wig.py -i {input} -o {bw_dir}/{wildcards.dataset} -t 1000000000 -s {chrom_size} -u  --strand='1+-,1-+,2++,2--'
        """

rule feature_count:
    input:
        expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam",work_dir = work_dir, dataset = samples)

    output:
        "{work_dir}/03_featurecount/F1_featureCounts.txt",
        "{work_dir}/03_featurecount/F1_featureStat.log"

    shell:
        "Rscript {script_dir}/featurecount.R"

rule sample_distance_plot:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)
 
    output:
        "{results_dir}/figure/distance_PNG_300_lengend_dpi.png",
        
    shell:
        "Rscript {script_dir}/distance.R"


rule DEgene_filter_and_plot_PSF:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/PSF/PSF_volcan_PNG_300_lengend_dpi.png",
        "{results_dir}/figure/PSF/DE_number.pdf",
        "{results_dir}/gene_name/PSF/up",
        "{results_dir}/gene_name/PSF/down",
        "{results_dir}/gene_name/PSF/sig",
        "{results_dir}/DEgene/PSF/ano_downgene.tsv",
        "{results_dir}/DEgene/PSF/ano_upgene.tsv"

    shell:
        "Rscript {script_dir}/PSF_DE.R"

rule DEgene_filter_and_plot_PSPC1:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/PSPC1/PSPC1_volcan_PNG_300_lengend_dpi.png",
        "{results_dir}/figure/PSPC1/DE_number.pdf",
        "{results_dir}/gene_name/PSPC1/up",
        "{results_dir}/gene_name/PSPC1/down",
        "{results_dir}/gene_name/PSPC1/sig",
        "{results_dir}/DEgene/PSPC1/ano_downgene.tsv",
        "{results_dir}/DEgene/PSPC1/ano_upgene.tsv"

    shell:
        "Rscript {script_dir}/PSPC1_DE.R"

rule DEgene_filter_and_plot_NONO:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/NONO/NONO_volcan_PNG_300_lengend_dpi.png",
        "{results_dir}/figure/NONO/DE_number.pdf",
        "{results_dir}/gene_name/NONO/up",
        "{results_dir}/gene_name/NONO/down",
        "{results_dir}/gene_name/NONO/sig",
        "{results_dir}/DEgene/NONO/ano_downgene.tsv",
        "{results_dir}/DEgene/NONO/ano_upgene.tsv"

    shell:
        "Rscript {script_dir}/NONO_DE.R"

rule DEgene_filter_and_plot_NEAT1_V2:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/NEAT1_V2/NEAT1_V2_volcan_PNG_300_lengend_dpi.png",
        "{results_dir}/figure/NEAT1_V2/DE_number.pdf",
        "{results_dir}/gene_name/NEAT1_V2/up",
        "{results_dir}/gene_name/NEAT1_V2/down",
        "{results_dir}/gene_name/NEAT1_V2/sig",
        "{results_dir}/DEgene/NEAT1_V2/ano_downgene.tsv",
        "{results_dir}/DEgene/NEAT1_V2/ano_upgene.tsv"

    shell:
        "Rscript {script_dir}/NEAT1_V2_DE.R"

rule DEgene_filter_and_plot_NEAT1_NT:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/NEAT1_NT/NEAT1_NT_volcan_PNG_300_lengend_dpi.png",
        "{results_dir}/figure/NEAT1_NT/DE_number.pdf",
        "{results_dir}/gene_name/NEAT1_NT/up",
        "{results_dir}/gene_name/NEAT1_NT/down",
        "{results_dir}/gene_name/NEAT1_NT/sig",
        "{results_dir}/DEgene/NEAT1_NT/ano_downgene.tsv",
        "{results_dir}/DEgene/NEAT1_NT/ano_upgene.tsv"

    shell:
        "Rscript {script_dir}/NEAT1_NT_DE.R"


