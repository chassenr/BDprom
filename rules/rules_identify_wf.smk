# specify config file
# change this path to the location of the config file with your project specific settings
configfile: "/srv/bio/Analysis_data/IOWseq000082_BDprom/Intermediate_results/Repos/BDprom/config/config_identify_wf.yaml"

rule all:
	input:
		ge_candidate = config["wdir"] + "/Pairs_out/bdprom_candidate_pairs.txt",
		ge_confirmed = config["wdir"] + "/Pairs_out/bdprom_ge_pairs.txt",
		ge_vis = config["wdir"] + "/Pairs_out/bdprom_plot_pairs.pdf"
	
rule identify_pairs:
	input:
		gtf = config["wdir"] + "/Genome/genome.gtf"
	output:
		ge_candidate = config["wdir"] + "/Pairs_out/bdprom_candidate_pairs.txt",
		ge_confirmed = config["wdir"] + "/Pairs_out/bdprom_ge_pairs.txt",
		ge_vis = config["wdir"] + "/Pairs_out/bdprom_plot_pairs.pdf"
	params:
		wdir = config["wdir"],
		dist_max = config["dist_max"],
		cor_min =  config["cor_min"],
		slope_diff = config["slope_diff"],
		median_diff = config["median_diff"],
		detect_f = config["detect_f"],
		sdir = config["sdir"]
	conda: "r-4.4.0"	
	log:
		config["wdir"] + "/Bdprom_pairs/logs/pairs.log"
	shell:
		"""
		mkdir -p {params.wdir}/Bdprom_pairs
		mkdir -p {params.wdir}/Bdprom_pairs/logs

		{params.sdir}/wrapper_identify_wf.R -s {params.sdir} -w {params.wdir} -d {params.dist_max} -c {params.cor_min} -l {params.slope_diff} -m {params.median_diff} -f {params.detect_f} -o {params.wdir}/Bdprom_pairs > {log} 2>&1
		"""
