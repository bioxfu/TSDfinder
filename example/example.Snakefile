configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R1_paired_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R2_paired_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}.bam', sample=config['samples']),
		expand('bam/{sample}.bam.bai', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
		expand('stat/bamqc_stat.tsv'),
		expand('track/{sample}.tdf', sample=config['samples']),
		expand('TE/index/TE_indexes_ok'),
		expand('output/{sample}', sample=config['samples']),

rule fastqc_raw_PE:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/raw {input}'

rule trimmomatic_PE:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_R1_paired.fastq.gz',
		'clean/{sample}_R2_paired.fastq.gz'
	output:
		'fastqc/clean/{sample}_R1_paired_fastqc.html',
		'fastqc/clean/{sample}_R2_paired_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R1_paired_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R2_paired_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule bowtie2_PE:
	input:
		r1 = 'clean/{sample}_R1_paired.fastq.gz',
		r2 = 'clean/{sample}_R2_paired.fastq.gz'
	output:
		bam = 'bam/{sample}.bam'
	params:
		prefix = 'bam/{sample}',
		cpu = config['cpu'],
		index = config['index'],
	shell:
		"bowtie2 -p {params.cpu} --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive -x {params.index} -1 {input.r1} -2 {input.r2} |samtools view -Shub|samtools sort - -T {params.prefix} -o {output.bam}"

rule bam_idx:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bai = 'bam/{sample}.bam.bai'
	shell:
		'samtools index {input.bam} {output.bai}'

rule bam_qc:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bamqc_dir = 'bam/{sample}.bamqc',
		bamqc_html = 'bam/{sample}.bamqc/qualimapReport.html'
	params:
		cpu = config['cpu']
	shell:
		"qualimap bamqc --java-mem-size=10G -nt {params.cpu} -bam {input.bam} -outdir {output.bamqc_dir}"

rule bam_qc_stat:
	input:
		['bam/{sample}.bamqc/qualimapReport.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/bamqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']		
	shell:
		"{params.Rscript} script/mapping_stat_by_bamqc.R"

rule bam2count:
	input:
		bam = 'bam/{sample}.bam'
	output:
		cnt = 'bam/{sample}.cnt'
	shell:
		"samtools view -c -F 4 {input.bam} > {output.cnt}"
		
rule bam2bedgraph:
	input:
		bam = 'bam/{sample}.bam',
		cnt = 'bam/{sample}.cnt'
	output:
		bg = 'track/{sample}.bedgraph'
	shell:
		"X=`awk '{{print 1/$1*1000000}}' {input.cnt}`; "
		"bedtools genomecov -ibam {input.bam} -bga -scale $X -split > {output.bg}"

rule bedgraph2tdf:
	input:
		bg = 'track/{sample}.bedgraph'
	output:
		tdf = 'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		"igvtools toTDF {input.bg} {output.tdf} {params.IGV}"

rule prepare_TE:
	input:
		TE_coord = config['TE_coord']
	output:
		'TE/TE_indexes_ok'
	params:
		genome = config['index']+'.fa'
	shell:
		'bash script/prepare_TEs.sh -i {input} -g {params.genome} -d TE/index'

rule SPLITREADER:
	input:
		bam = 'bam/{sample}.bam',
		TE_coord = config['TE_coord']
	output:
		'output/{sample}'
	params:
		l = config['max_len'],
		r = config['min_split_reads'],
		m = config['maxcov'],
		p = config['cpu'],
		g = config['index']
	shell:
		'bash script/SPLITREADER.sh -o {output} -t {output} -b {input.bam} -e {input.TE_coord} -l {params.l} -r {params.r} -m {params.m} -p {params.p} -s TE/index -g {params.g}'