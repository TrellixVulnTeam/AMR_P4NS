assembly_K_list = "27,47,67,87,107,127"
rgi_annotation_alignment_tool = "BLAST"
k_mer_size = 27


		
rule all:
	input:
		"{sample}/read_classification.tsv",
		"{sample}/kmc_out/kmer_count_1.txt",
		"{sample}/kmc_out/kmer_count_2.txt",
		"{sample}/rgi_annotation.json",
		"{sample}/bwt.allele_mapping_data.json"
	output:
		"{sample}/analisi_completa.log"
	conda:"env.yaml"
	shell:
		"""
		touch {wildcards.sample}/analisi_completa.log
		"""		

rule population_set_up:
	input:
		i_1 = "{sample}/input.json",
		i_2 = "{sample}/input.tsv",
		g = "{sample}/genomes"
	output:
		"{sample}/config.ini"
	shell:
		"""
		python my_scripts/generate_files.py {wildcards.sample}
		"""

rule camisim:
	input:
		C = "{sample}/config.ini"
	output:
		"{sample}/reads_1.fastq",
		"{sample}/reads_2.fastq"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		rm -rf {wildcards.sample}/camisim_out
		mkdir -p {wildcards.sample}/tmp
		python CAMISIM/metagenomesimulation.py {input.C}
		gunzip {wildcards.sample}/camisim_out/2022*/reads/anonymous_reads.fq.gz
		python my_scripts/split_pair_end.py {wildcards.sample}/camisim_out/2022*/reads/anonymous_reads.fq {wildcards.sample}/
		"""

rule assembly:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq"
	output:
		"{sample}/assembly/spades/contigs.fasta",
		"{sample}/assembly/assembly_report.tsv"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		mkdir -p {wildcards.sample}/assembly
		spades.py -1 {input.S1} -2 {input.S2} -t {threads} -m 200 -k {assembly_K_list} -o {wildcards.sample}/assembly/spades --meta
		python my_scripts/assembly_report.py {wildcards.sample}/assembly/spades/contigs.fasta {wildcards.sample}/assembly/assembly_report.tsv
		"""

rule taxonomy:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq",
		db = "kraken2_stuff/KRAK_DB"
	output:
		"{sample}/taxonomy/read_classification.tsv",
		"{sample}/taxonomy/classified_reads_1.fastq",
		"{sample}/taxonomy/classified_reads_2.fastq",
		"{sample}/taxonomy/read_species.tsv",
		"{sample}/taxonomy/read_taxonomy.tsv"
	conda: "env.yaml"
	threads: 4
	shell:
		"""
		mkdir -p {wildcards.sample}/taxonomy
		kraken2 --db {input.db} --threads {threads} --classified-out {wildcards.sample}/taxonomy/classified_reads_#.fastq --paired {input.S1} {input.S2} > {wildcards.sample}/taxonomy/read_classification.tsv
		kraken-translate {wildcards.sample}/taxonomy/read_classification.tsv {wildcards.sample}/taxonomy/read_species.tsv
		kraken-report {wildcards.sample}/taxonomy/read_classification.tsv {wildcards.sample}/taxonomy/read_taxonomy.tsv
		"""


rule kmer_count:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq"
	output:
		"{sample}/kmc_out/kmer_count_1.txt.gz",
		"{sample}/kmc_out/kmer_count_2.txt.gz"
	conda: "env.yaml"
	threads: 20
	params:
		ram = 10
	shell:
		"""
		rm -r {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out
		mkdir {wildcards.sample}/kmc_out/kmc_temp
		kmc -k{k_mer_size} -m{params.ram} -t{threads} {input.S1} {wildcards.sample}/kmc_out/k_1 {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/k_1 dump {wildcards.sample}/kmc_out/kmer_count_1.txt
		gzip {wildcards.sample}/kmc_out/kmer_count_1.txt
		kmc -k{k_mer_size} -m{params.ram} -t{threads} {input.S2} {wildcards.sample}/kmc_out/k_2 {wildcards.sample}/kmc_out/kmc_temp
		kmc_tools transform {wildcards.sample}/kmc_out/k_2 dump {wildcards.sample}/kmc_out/kmer_count_2.txt
		gzip {wildcards.sample}/kmc_out/kmer_count_2.txt
		"""



rule rgi_annotation:
	input:
		C = "{sample}/assembly/spades/contigs.fasta",
		r_1 = "rgi/card_annotation.log",
		r_2 = "rgi/wildcard_annotation.log"
	output:
		"{sample}/rgi_annotation.json"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		rgi main -i {input.C} -o {wildcards.sample}/rgi_annotation -t contig -a {rgi_annotation_alignment_tool} --clean -n {threads} -d wgs --split_prodigal_jobs
		"""

rule rgi_bowtie2:
	input:
		S1 = "{sample}/reads_1.fastq",
		S2 = "{sample}/reads_2.fastq",
		r1 = "rgi/card_annotation.log",
		r2 = "rgi/wildcard_annotation.log"
	output:
		"{sample}/bwt_out/bwt.allele_mapping_data.json"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		mkdir -p {wildcards.sample}/bwt_out
		cd rgi
		rgi bwt -1 ../{input.S1} -2 ../{input.S2} -a bowtie2 -n {threads} --clean --include_wildcard -o ../{wildcards.sample}/bwt_out/bwt --local 
		cd ..
		"""

rule prepair_data:
	input:
		bwt = "{sample}/bwt_out/bwt.allele_mapping_data.json",
		ab = "{sample}/abundances.tsv",
		met = "{sample}/metadata.json"
	output:
		"{sample}/to_download.zip"
	shell:
		"""
		mkdir -p {wildcards.sample}/to_download
		cp {input.bwt} {input.ab} {input.meta} {wildcards.sample}/to_download/
		zip -r to_download.zip to_download
		rm -r to_download
		"""

rule kraken_env_setup:
	output:
		"kraken2_stuff/KRAK_DB"
	conda: "env.yaml"
	threads: 20
	shell:
		"""
		mkdir -p kraken_stuff
		kraken2-build --standard --threads {threads} --db kraken_stuff/KRAK_DB
		"""

rule rgi_env_setup:
	output:
		"rgi/card_annotation.log",
		"rgi/wildcard_annotation.log"
	conda: "env.yaml"
	shell:
		"""
		mkdir -p rgi
		cd rgi
		#CARD DATABASE
		wget -O data https://card.mcmaster.ca/latest/data
		tar -xvf data ./card.json
		rgi load --card_json card.json --local
		rgi card_annotation -i card.json > card_annotation.log 2>&1
		rgi load -i card.json --card_annotation card_database* --local
		#WILDCARD DATABASE
		wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
		mkdir -p wildcard
		tar -xjf wildcard_data.tar.bz2 -C wildcard
		gunzip wildcard/*.gz
		rgi wildcard_annotation -i wildcard --card_json card.json -v_non_specificata > wildcard_annotation.log 2>&1
		rgi load --wildcard_annotation wildcard_database* --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database* --local
		cd ..
		"""
