# specify config file
# change this path to the location of the config file with your project specific settings
configfile: "/otc-genomics/IOWseq000082_BDprom/Intermediate_results/Repos/tmp/BDprom/config/config_identify_wf.yaml"

rule all:
	input:
		ge_pairs = config["wdir"] + "/Bdprom_pairs/bdprom_candidate_pairs.txt"
	
rule identifygenepairs:
	input:
		gtf = config["wdir"] + "/Genome/genome.gtf"
	output:
		ge_pairs = config["wdir"] + "/Bdprom_pairs/bdprom_candidate_pairs.txt"
	params:
		script = "src/identify_confirm_bdgenepairs.R"
	conda: "r-4.4.0"	
	log:
		config["wdir"] + "/Bdprom_pairs/logs/pairs.log"
	shell:
		"""
		{params.script} {output.ge_pairs} > {log} 2>&1	
		"""
