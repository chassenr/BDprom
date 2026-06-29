# specify config file
# change this path to the location of the config file with your project specific settings
configfile: "/otc-genomics/IOWseq000082_BDprom/Intermediate_results/Repos/BDprom/config/config_identify_wf.yaml"

rule all:
	input:
		ge_pairs = config["wdir"] + "/Bdprom_pairs/bdprom_candidate_pairs.txt"
	
rule identifygenepairs:
	input:
		gtf = config["wdir"] + "/Phaeodactylum_tricornutum/Genome/genome.gtf"
	output:
		ge_pairs = config["wdir"] + "/Bdprom_pairs/bdprom_candidate_pairs.txt"
	params:
		indir =  config["wdir"] + "/Phaeodactylum_tricornutum/Parsed_counts",
		cor_min =  config["cor_min"],
		ge_dir = config["ge_dir"],
		slope_diff = config["slope_diff"],
		median_diff = config["median_diff"],
		detect_f = config["detect_f"],
		script = "Intermediate_results/identify_confirm_bdgenepairs.R"
	conda: "r-4.4.0"	
	log:
		config["wdir"] + "/Bdprom_pairs/logs/pairs.log"
	shell:
		"""
		-s {params.script} -i {params.indir} -c {params.indir} -g {params.ge_dir} -l {params.slope_diff} -m  {params.median_diff} -t {params.detect_f}  {input.gtf} -o {output.ge_pairs} > {log} 2>&1	
		"""
