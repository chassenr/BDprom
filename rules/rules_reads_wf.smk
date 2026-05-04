# specify config file
# change this path to the location of the config file with your project specific settings
configfile: "/otc-genomics/IOWseq000082_BDprom/Intermediate_results/Repos/BDprom/config/config_reads_wf.yaml"

rule all:
    input: 
        ena_out = config["wdir"] + "/ENA_out/RNAseq_raw.tsv",
        ena_filtered = config["wdir"] + "/ENA_out/RNAseq_filtered.tsv"

rule search_ena:
    output:
        ena_out = config["wdir"] + "/ENA_out/RNAseq_raw.tsv"
    params:
        taxid = config["tax_id"]
    shell:
        """
        curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_run&query=tax_tree('{params.tax_id}')%20AND%20library_source%3D%22transcriptomic%22%20AND%20library_strategy%3D%22rna-seq%22&fields=run_accession%2Cexperiment_title%2Ctax_id%2Clibrary_source%2Cbase_count%2Cbio_material%2Cbroad_scale_environmental_context%2Cbroker_name%2Ccenter_name%2Cchecklist%2Ccollected_by%2Ccollection_date%2Ccollection_date_end%2Ccollection_date_start%2Ccountry%2Cdescription%2Cenvironmental_medium%2Cexperiment_accession%2Cexperiment_alias%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_md5%2Cfirst_created%2Cfirst_public%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_construction_protocol%2Clibrary_gen_protocol%2Clibrary_layout%2Clibrary_max_fragment_size%2Clibrary_min_fragment_size%2Clibrary_name%2Clibrary_selection%2Clibrary_strategy%2Clocal_environmental_context%2Clocation%2Clocation_end%2Clocation_start%2Clon%2Cmarine_region%2Cncbi_reporting_standard%2Cnominal_length%2Cnominal_sdev%2Cproject_name%2Cread_count%2Cread_strand%2Crna_integrity_num%2Crna_purity_230_ratio%2Crna_purity_280_ratio%2Crun_alias%2Csample_accession%2Csample_alias%2Csample_collection%2Csample_description%2Csample_material%2Csample_storage%2Csample_storage_processing%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_project%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_method%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmission_tool%2Ctag%2Ctarget_gene%2Ctax_lineage%2Ctaxonomic_classification&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search" > {output.ena.out}
        """

rule filter_metadata:
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
        config["wdir"] + "ENA_out/logs/filter_ena_metadata.log"
    shell:
        """
        {params.script} -e {input.ena_curl} -c {params.min_read_count} -r {params.min_n_runs} -s {params.seq_technology} -p {params.only_keep_pubmed} -o {output.ena_filtered} >> {log} 2>&1
        """
    
    rule download_files:
        input:
        output:
        log:
        shell:
            """

            """