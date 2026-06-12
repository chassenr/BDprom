# To import the RNAseq_filteredR config["wdir"] + "/Parsed_counts/NAseq_filtered.tsv file two packages are needed:
import glob
import pandas as pd

# specify config file
# change this path to the location of the config file with your project specific settings
configfile: "/otc-genomics/IOWseq000082_BDprom/Intermediate_results/Repos/tmp/BDprom/config/config_reads_wf.yaml"

# def aggregate_study (wildcards):
#	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
#	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n")
#	STUDY = RNAseq_filtered.study_accession.unique().to_list()
#	return expand(config["wdir"] + "/Parsed_counts/gene_expression_{study}.csv", study = STUDY)

rule all:
	input:
		ena_out = config["wdir"] + "/ENA_out/RNAseq_raw.tsv",
		ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv",
		pmc = config["wdir"] + "/Abstracts/pmc_abstracts.txt" if config["only_keep_pubmed"] == "TRUE" else [],
		pubmed = config["wdir"] + "/Abstracts/pubmed_abstracts.txt" if config["only_keep_pubmed"] == "TRUE" else [],
		done = config["wdir"] + "/Fastq_downloaded/done",
		fastqc = config["wdir"] + "/Fastq_qc/done",
		parsed_counts = config["wdir"] + "/Parsed_counts/done",
		summary_phix = config["wdir"] + "/Logs_summary/summary_phix.tsv",
		summary_counts = config["wdir"] + "/Logs_summary/summary_counts.csv"

rule search_ena:
	output:
		ena_out = config["wdir"] + "/ENA_out/RNAseq_raw.tsv"
	params:
		tax_id = config["tax_id"]
	shell:
		"""
		curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree({params.tax_id})%20AND%20library_source%3D%22transcriptomic%22%20AND%20library_strategy%3D%22rna-seq%22&fields=run_accession%2Cexperiment_title%2Ctax_id%2Clibrary_source%2Cbase_count%2Cbio_material%2Cbroad_scale_environmental_context%2Cbroker_name%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccollection_date_end%2Ccollection_date_start%2Ccountry%2Cdescription%2Cenvironmental_medium%2Cexperiment_accession%2Cexperiment_alias%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_construction_protocol%2Clibrary_gen_protocol%2Clibrary_layout%2Clibrary_max_fragment_size%2Clibrary_min_fragment_size%2Clibrary_name%2Clibrary_selection%2Clibrary_strategy%2Clocal_environmental_context%2Clocation%2Clocation_end%2Clocation_start%2Clon%2Cmarine_region%2Cncbi_reporting_standard%2Cnominal_length%2Cnominal_sdev%2Cproject_name%2Cread_count%2Cread_strand%2Crna_integrity_num%2Crna_purity_230_ratio%2Crna_purity_280_ratio%2Crun_alias%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_storage%2Csample_storage_processing%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_project%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmission_tool%2Ctag%2Ctarget_gene%2Ctax_lineage%2Ctaxonomic_classification&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" > {output.ena_out}
		"""

checkpoint filter_metadata:
	input:
		ena_out = config["wdir"] + "/ENA_out/RNAseq_raw.tsv"
	output:
		ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"
	params:
		min_read_counts = config["min_read_counts"],
		min_n_runs = config["min_n_runs"],
		seq_technology = config["seq_technology"],
		only_keep_pubmed = config["only_keep_pubmed"],
		script = "src/filter_ena_metadata.R"
	conda: "r-4.4.0"
	log:
		config["wdir"] + "/ENA_out/logs/filter_ena_metadata.log"
	shell:
		"""
		{params.script} -e {input.ena_out} -c {params.min_read_counts} -r {params.min_n_runs} -s {params.seq_technology} -p {params.only_keep_pubmed} -o {output.ena_filtered} >> {log} 2>&1
		"""

if config["only_keep_pubmed"] == "TRUE":
	rule get_abstract_pmc:
		input: 
			ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"
		output:
			pmc = config["wdir"] + "/Abstracts/pmc_abstracts.txt"
		params:
			outdir = config["wdir"] + "/Abstracts/tmp_pmc",
			script = "src/parse_abstracts.sh"
		conda: "anvio-8"
		shell:
			"""
			mkdir -p {params.outdir}
			echo -E "pmc_id\ttitle\tabstract" > {output.pmc}
			CN=$(head -1 {input.ena_filtered} | tr '\\t' '\\n' | grep -i -w -n "EuropePMC" | cut -d':' -f1)
			cut -f$CN {input.ena_filtered} | sed '1d' | tr ';' '\n' | sort -u | sed '/^NA$/d' | while read line
			do
		  	  efetch -db pmc -format medline -mode text -id $line > {params.outdir}/$line.efetch
			  {params.script} {params.outdir}/$line.efetch >> {output.pmc}
			done
			"""

if config["only_keep_pubmed"] == "TRUE":
	rule get_abstract_pubmed:
		input: 
			ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"
		output:
			pubmed = config["wdir"] + "/Abstracts/pubmed_abstracts.txt"
		params:
			outdir = config["wdir"] + "/Abstracts/tmp_pubmed",
			script = "src/parse_abstracts.sh"
		conda: "anvio-8"
		shell:
			"""
			mkdir -p {params.outdir}
			echo -E "pubmed_id\ttitle\tabstract" > {output.pubmed}
			CN=$(head -1 {input.ena_filtered} | tr '\\t' '\\n' | grep -i -w -n "PubMed" | cut -d':' -f1)
			cut -f$CN {input.ena_filtered} | sed '1d' | tr ';' '\n' | sort -u | sed '/^NA$/d' | while read line
			do
		  	  efetch -db pubmed -format medline -mode text -id $line > {params.outdir}/$line.efetch
		  	  {params.script} {params.outdir}/$line.efetch >> {output.pubmed}
			done
			"""

rule download_files:
	input:
		ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"
	output:
		done = config["wdir"] + "/Fastq_downloaded/done"
	params:
		outdir = config["wdir"] + "/Fastq_downloaded"
	conda: "foldseek-10.941cd33"
	log:
		aria2 =config["wdir"] + "/Fastq_downloaded/logs/aria2c.log",
		stdout = config["wdir"] + "/Fastq_downloaded/logs/md5_out.log",
		stderr = config["wdir"] + "/Fastq_downloaded/logs/md5_err.log"
	shell:
		"""
		set -x
		CN=$(head -1 {input.ena_filtered} | tr '\\t' '\\n' | grep -w -n "fastq_ftp" | cut -d':' -f1)
		cut -f$CN {input.ena_filtered} | sed '1d' | tr ';' '\\n' | sed 's/^/http:\\/\\//' > {params.outdir}/download_links.txt
		(aria2c -i {params.outdir}/download_links.txt -c --dir {params.outdir} --max-tries=20 --retry-wait=5 --max-connection-per-server=1 --max-concurrent-downloads=16 -V >> {log.aria2} 2>&1) || true
		MD5=$(head -1 {input.ena_filtered} | tr '\\t' '\\n' | grep -w -n "fastq_md5" | cut -d':' -f1)
		paste <(cut -f$MD5 {input.ena_filtered} | sed '1d' | tr ';' '\\n') <(cut -f$CN {input.ena_filtered} | sed '1d' | tr ';' '\\n' | xargs -n1 basename | sed "s|^|{params.outdir}/|") > {params.outdir}/md5_checksums.txt
		(md5sum -c {params.outdir}/md5_checksums.txt > {log.stdout} 2> {log.stderr}) || true
		if [ $(wc -l < {params.outdir}/md5_checksums.txt) == $(grep -c "OK" {log.stdout}) ]
		then
		  echo "All downloaded files passed checksums."
		  touch {output.done}
		else
		  echo "Not all files were downloaded correctly. Please remove them and repeat the download.\n"
		  grep "FAILED" {log.stdout}
		fi
		"""

rule phix_clipping:
	input:
		done = config["wdir"] + "/Fastq_downloaded/done"
	output:
		phix_clipped = config["wdir"] + "/Fastq_phix_clipped/{sample}_1.fastq.gz"
	params:
		indir = config["wdir"] + "/Fastq_downloaded",
		outdir = config["wdir"] + "/Fastq_phix_clipped",
		sample = "{sample}"
	threads: config["parallel_threads"]
	conda: "busco-5.8.0"
	log:
		config["wdir"] + "/Fastq_phix_clipped/logs/{sample}_phix.log"
	shell:
		"""
		if [ -f {params.indir}/{params.sample}_2.fastq.gz ]
		then
		  bbduk.sh in={params.indir}/{params.sample}_1.fastq.gz in2={params.indir}/{params.sample}_2.fastq.gz out={output.phix_clipped} out2={params.outdir}/{params.sample}_2.fastq.gz ref=phix k=28 threads={threads} > {log} 2>&1
		else
		  bbduk.sh in={params.indir}/{params.sample}.fastq.gz out={output.phix_clipped} ref=phix k=28 threads={threads} > {log} 2>&1
		fi
		"""

rule adapter_clipping:
	input:
		phix_clipped = config["wdir"] + "/Fastq_phix_clipped/{sample}_1.fastq.gz"
	output:
		adapter_clipped = config["wdir"] + "/Fastq_adapter_clipped/{sample}_1.fastq.gz"
	params:
		indir = config["wdir"] + "/Fastq_phix_clipped",
		outdir = config["wdir"] + "/Fastq_adapter_clipped",
		sample = "{sample}"
	threads: config["parallel_threads"]
	conda: "busco-5.8.0"
	log:
		config["wdir"] + "/Fastq_adapter_clipped/logs/{sample}_adapter.log"
	shell:
		"""
		if [ -f {params.indir}/{params.sample}_2.fastq.gz ]
		then
		  bbduk.sh in={input.phix_clipped} in2={params.indir}/{params.sample}_2.fastq.gz out={output.adapter_clipped} out2={params.outdir}/{params.sample}_2.fastq.gz ref=adapters k=23 mink=11 ktrim=r hdist=1 threads={threads} > {log} 2>&1
		else
		  bbduk.sh in={input.phix_clipped} out={output.adapter_clipped} ref=adapters k=23 mink=11 ktrim=r hdist=1 threads={threads} > {log} 2>&1  
		fi
		"""
rule trimming:
	input:
		adapter_clipped = config["wdir"] + "/Fastq_adapter_clipped/{sample}_1.fastq.gz"
	output:
		trimmed = config["wdir"] + "/Fastq_trimmed/{sample}_1.fastq.gz"
	params:
		indir = config["wdir"] + "/Fastq_adapter_clipped",
		outdir = config["wdir"] + "/Fastq_trimmed",
		sample = "{sample}"
	threads: config["parallel_threads"]
	conda: "busco-5.8.0" 
	log:
		trimmed = config["wdir"] + "/Fastq_trimmed/logs/{sample}_trimmed.log"
	shell:
		"""
		# set -x
		MAXLEN=$(zcat {input.adapter_clipped} | head -4000 | sed -n '2~4p' | perl -nle 'print length' | sort -k1,1gr | head -1 || true)
		MINLEN=$(bc <<< "${{MAXLEN}} * 3 / 4")
		if [ -f {params.indir}/{params.sample}_2.fastq.gz ]
		then
		  bbduk.sh in={input.adapter_clipped} in2={params.indir}/{params.sample}_2.fastq.gz out={output.trimmed} out2={params.outdir}/{params.sample}_2.fastq.gz qtrim=w,4 trimq=15 minlength=$MINLEN trimpolyg=10 threads={threads} > {log} 2>&1
		else
		  bbduk.sh in={input.adapter_clipped} out={output.trimmed} qtrim=w,4 trimq=15 minlength=$MINLEN trimpolyg=10 threads={threads} > {log} 2>&1
		fi
		"""

def aggregate_fastqc (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Fastq_trimmed/{sample}_1.fastq.gz", sample = SAMPLES)
		
rule fastqc:
	input:
		trimmed = aggregate_fastqc
	output: 
		done = config["wdir"] + "/Fastq_qc/done"
	params:
		indir = config["wdir"] + "/Fastq_trimmed",
		outdir = config["wdir"] + "/Fastq_qc"
	threads: config["parallel_threads"]
	conda: "qc-env"
	shell:
		"""
		fastqc -t {threads} -o {params.outdir} {params.indir}/*.fastq.gz
		touch {output.done}
		"""

rule genome_download:
	output:
		gtf = config["wdir"] + "/Genome/genome.gtf",
		genome = config["wdir"] + "/Genome/genome.fa"
	params:
		outdir = config["wdir"] + "/Genome",
		gtf = config["gtf_link"],
		fasta = config["genome_link"],
		tmp_fa = config["wdir"] + "/Genome/genome.fa.gz",
		tmp_gtf = config["wdir"] + "/Genome/genome.gtf.gz"
	conda: "foldseek-10.941cd33"
	shell:
 		"""
		wget -O {params.tmp_gtf} {params.gtf}
		wget -O {params.tmp_fa} {params.fasta}
		gunzip {params.tmp_gtf} {params.tmp_fa} 
		"""

rule index:
	input:
		genome = config["wdir"] + "/Genome/genome.fa",
		gtf = config["wdir"] + "/Genome/genome.gtf"
	output:
		splice = config["wdir"] + "/Genome/index/genome_ss.txt",
		exon = config["wdir"] + "/Genome/index/genome_exon.txt",
		index = config["wdir"] + "/Genome/index/genome.1.ht2"
	params:
		index = config["wdir"] + "/Genome/index/genome"
	conda: "braker-3.0.8"
	log:
		config["wdir"] + "/Genome/index/genome.log"
	shell:
		"""
		hisat2_extract_exons.py {input.gtf} > {output.exon}
		hisat2_extract_splice_sites.py {input.gtf} > {output.splice}
		hisat2-build --exon {output.exon} --ss {output.splice} -f {input.genome} {params.index} > {log} 2>&1
		"""

rule mapping:
	input:
		trimmed = config["wdir"] + "/Fastq_trimmed/{sample}_1.fastq.gz",
		index = config["wdir"] + "/Genome/index/genome.1.ht2"
	output:
		bam = config["wdir"] + "/Fastq_mapped/{sample}.bam"
	params:
		indir = config["wdir"] + "/Fastq_trimmed",
		sample = "{sample}",
		outdir = config["wdir"] + "/Fastq_mapped",
		index = config["wdir"] + "/Genome/index/genome"
	conda: "braker-3.0.8"
	threads: config["parallel_threads"]
	log:
		config["wdir"] + "/Fastq_mapped/logs/{sample}_mapped.log"
	shell:
		"""
		if [ -f {params.indir}/{params.sample}_2.fastq.gz ]
		then
		  hisat2 -p {threads} --no-unal -x {params.index} -1 {input.trimmed} -2 {params.indir}/{params.sample}_2.fastq.gz 2>{log} | samtools sort -T {params.outdir}/tmp_{params.sample} -@{threads} | samtools view -F2304 -f2 -b -@{threads} -o {output.bam}
		else
		  hisat2 -p {threads} --no-unal -x {params.index} -U {input.trimmed} 2>{log} | samtools sort -T {params.outdir}/tmp_{params.sample} -@{threads} | samtools view -F2304 -b -@{threads} -o {output.bam}
		fi
		"""

rule feature_counts:
	input:
		bam = config["wdir"] + "/Fastq_mapped/{sample}.bam",
		gtf = config["wdir"] + "/Genome/genome.gtf"
	output:
		counts = config["wdir"] + "/Gene_expression/{sample}_counts.txt"
	params:
		trimdir = config["wdir"] + "/Fastq_trimmed",
		sample = "{sample}"		
	conda: "subread-2.0.6"
	threads: config["parallel_threads"] 
	log:
		config["wdir"] + "/Gene_expression/logs/{sample}_counts.log"
	shell:
		"""
		if [ -f {params.trimdir}/{params.sample}_2.fastq.gz ]
		then
		  featureCounts -t exon -g gene_id -p --countReadPairs --fracOverlap 0.8 -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		else
		  featureCounts -t exon -g gene_id --fracOverlap 0.8 -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		fi
		"""

def aggregate_counts (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Gene_expression/{sample}_counts.txt", sample = SAMPLES)

rule parse_counts:
	input:
		counts = aggregate_counts,
		ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"
	output:
		done = config["wdir"] + "/Parsed_counts/done"
	params:
		indir = config["wdir"] + "/Gene_expression",
		outdir = config["wdir"] + "/Parsed_counts",
		script = "src/parse_feature_counts.R"
	conda: "r-4.4.0"
	log:
		config["wdir"] + "/Parsed_counts/logs/parsed.log"
	shell:
		"""
		{params.script} -i {params.indir} -o {params.outdir} -e {input.ena_filtered} >> {log} 2>&1	
		touch {output.done}
		"""

def aggregate_phix (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Fastq_phix_clipped/logs/{sample}_phix.log", sample = SAMPLES)

def aggregate_adapter (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Fastq_adapter_clipped/logs/{sample}_adapter.log", sample = SAMPLES)

def aggregate_trimmed (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Fastq_trimmed/logs/{sample}_trimmed.log", sample = SAMPLES)

def aggregate_mapped (wildcards):
	Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
	RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n") 
	SAMPLES = RNAseq_filtered.run_accession.to_list()
	return expand(config["wdir"] + "/Fastq_mapped/logs/{sample}_mapped.log", sample = SAMPLES)

def aggregate_counts (wildcards):
        Check_rnaseq_filtered = checkpoints.filter_metadata.get().output.ena_filtered
        RNAseq_filtered = pd.read_table(Check_rnaseq_filtered, sep = "\t", header = 0, lineterminator = "\n")
        SAMPLES = RNAseq_filtered.run_accession.to_list()
        return expand(config["wdir"] + "/Gene_expression/logs/{sample}_counts.log", sample = SAMPLES)

rule collect_logs:
	input:
		phix = aggregate_phix,
		adapter = aggregate_adapter, 
		trimmed = aggregate_trimmed, 
		mapped = aggregate_mapped,
		counts = aggregate_counts
	output:
		summary_phix = config["wdir"] + "/Logs_summary/summary_phix.tsv",
		summary_adapter = config["wdir"] + "/Logs_summary/summary_adapter.tsv",
		summary_trim = config["wdir"] + "/Logs_summary/summary_trimmed.tsv",
		summary_mapped = config["wdir"] + "/Logs_summary/summary_mapped.tsv",
		summary_counts = config["wdir"] + "/Logs_summary/summary_counts.csv"
	shell:
		"""
		echo -e 'run_accession\\treads\\treads_perc\\tbp\\tbp_perc' > {output.summary_phix} 
		grep "Result" {input.phix} | sed -E 's/\\s+/\\t/g' | sed -e 's/_phix\\.log:Result://' -e 's/^.*\\///' | cut -f1,2,4,5,7 | sed -e 's/(//g' -e 's/%)//g' >> {output.summary_phix}
		echo -e 'run_accession\\treads\\treads_perc\\tbp\\tbp_perc' > {output.summary_adapter}
		grep "Result" {input.adapter} | sed -E 's/\\s+/\\t/g' | sed -e 's/_adapter\\.log:Result://g' -e 's/^.*\\///' | cut -f1,2,4,5,7 | sed -e 's/(//g' -e 's/%)//g' >> {output.summary_adapter}
		echo -e 'run_accession\\treads\\treads_perc\\tbp\\tbp_perc' > {output.summary_trim}
		grep "Result" {input.trimmed} | sed -E 's/\\s+/\\t/g' | sed -e 's/_trimmed\\.log:Result://g' -e 's/^.*\\///' | cut -f1,2,4,5,7 | sed -e 's/(//g' -e 's/%)//g' >> {output.summary_trim}
		echo -e 'run_acession\\toverall_alignment_perc' > {output.summary_mapped}
		grep "overall alignment rate" {input.mapped} | sed -e 's/_mapped\\.log:/ /' -e 's/%//' -e 's/^.*\\///' | sed -E 's/\\s+/\\t/g' | cut -f1,2 >> {output.summary_mapped}
		
		echo -e 'run_accession\treads\talignment_perc' > {output.summary_counts}
		grep "Successfully assigned alignments" {input.counts} | sed -E 's/\\s+/\\t/g' | sed -e 's/_counts.log:||//' -e 's/(//g' -e 's/%)//' -e 's/^.*\\///' | cut -f1,6,7 >> {output.summary_counts}
		"""
