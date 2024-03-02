import groovy.json.JsonOutput
/*
* This process takes the output produced by the
* `salmon_map_rad` rule and generates a permit list
* with alevin-fry (currently using the knee method).
*/
process alevin_fry_gpl {
    tag "fry_gpl:${fastq_MD5sum}:${filt_type}"
    // conda "bioconda::alevin-fry"
    label 'single_threads'

    input:
        val chemistry
        val reference
        val dataset_name
        val dataset_url
        val fastq_url
        val fastq_MD5sum
        val delete_fastq
        val feature_barcode_csv_url
        val multiplexing_library_csv_url
        val filt_type
        path pl_path
        path t2g_path
        path map_dir_path

    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_url, emit:fastq_url
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        val filt_type, emit: filt_type
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path
        path map_dir_path, emit: map_dir_path
        path "${fastq_MD5sum}_permitlist_${filt_type}_fw", emit: permit_dir_path

    script:
        filt_flag = (filt_type == "knee") ? "-k" : "-u ${pl_path}"

        cmd = """
        $params.timecmd -v -o ${fastq_MD5sum}_logs_permitlist_${filt_type}.time alevin-fry generate-permit-list ${filt_flag} -d fw -i ${map_dir_path} -o ${fastq_MD5sum}_permitlist_${filt_type}_fw
        """

        """
        ${cmd}
        """

    stub:
        filt_flag = (filt_type == "knee") ? "-k" : "-u ${pl_path}"
        cmd = """
        $params.timecmd -v -o ${fastq_MD5sum}_logs_permitlist_${filt_type}.time alevin-fry generate-permit-list ${filt_flag} -d fw -i ${map_dir_path} -o ${fastq_MD5sum}_permitlist_${filt_type}_fw
        """

        """
        echo "running ${cmd}"
        echo ${map_dir_path}
        echo ${fastq_MD5sum}_permitlist_${filt_type}_fw
        mkdir ${fastq_MD5sum}_permitlist_${filt_type}_fw
        touch ${fastq_MD5sum}_permitlist_${filt_type}_fw/permit_map.bin
        """
}

/*
* This process takes the output produced by the
* `alevin_fry_gpl` rule and generates a collated
* RAD file.
*/
process alevin_fry_collate {
    tag "fry_collate:${fastq_MD5sum}:${filt_type}"
    label 'multi_threads'
    // conda "bioconda::alevin-fry"

    input:
        val chemistry
        val reference
        val dataset_name
        val dataset_url
        val fastq_url
        val fastq_MD5sum
        val delete_fastq
        val feature_barcode_csv_url
        val multiplexing_library_csv_url
        val filt_type
        path pl_path
        path t2g_path
        path map_dir_path
        path permit_dir_path

    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_url, emit:fastq_url
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        val filt_type, emit: filt_type
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path
        path map_dir_path, emit: map_dir_path
        path permit_dir_path, emit: collate_dir

    script:
        cmd = """
        $params.timecmd -v -o ${fastq_MD5sum}_logs_collate_${filt_type}.time alevin-fry collate \
    -i ${permit_dir_path} -r ${map_dir_path} -t ${task.cpus}
        """

        """
        ${cmd}
        """

    stub:
    cmd = """
        $params.timecmd -v -o ${fastq_MD5sum}_logs_collate_${filt_type}.time alevin-fry collate \
    -i ${permit_dir_path} -r ${map_dir_path} -t ${task.cpus}
    """

    """
    echo "executing :: ${cmd}"
    # mkdir -p ${permit_dir_path}
    touch ${permit_dir_path}/map.collated.rad
    """
}

/*
* This process takes the output produced by the
* `alevin_fry_collate` rule and generates the resulting
* quantification.
*/
process alevin_fry_quant {
    tag "fry_quant:${fastq_MD5sum}:${filt_type}"
    // conda "bioconda::alevin-fry"
    label 'multi_threads'

    publishDir "${params.output_dir}/alevin_fry", mode: 'copy'

    input:
        val chemistry
        val reference
        val dataset_name
        val dataset_url
        val fastq_url
        val fastq_MD5sum
        val delete_fastq
        val feature_barcode_csv_url
        val multiplexing_library_csv_url
        val filt_type
        path pl_path
        path t2g_path
        path map_dir_path
        path collate_dir

    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_url, emit:fastq_url
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        val filt_type, emit: filt_type
        path "${odir}/**", emit: quant_dir

    script:
        odir = "${fastq_MD5sum}_fry_${filt_type}_quant_usa_cr-like"
        cmd = """$params.timecmd -v -o ${fastq_MD5sum}_logs_quant_${filt_type}.time alevin-fry quant \
    -r cr-like --use-mtx -m ${t2g_path} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""
        """
            ${cmd}
        """

    stub:
        odir = "${fastq_MD5sum}_fry_${filt_type}_quant_usa_cr-like"
        cmd = """$params.timecmd -v -o ${fastq_MD5sum}_logs_quant_${filt_type}.time alevin-fry quant \
    -r cr-like --use-mtx -m ${t2g_path} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""

        """
            echo "executing :: ${cmd}"
            mkdir -p ${odir}
            touch ${odir}/meta_info.json
        """
}

process write_description {
    tag "write_description:${fastq_MD5sum}:${filt_type}"
    // conda "bioconda::alevin-fry"
    label 'single_threads'

    publishDir "${params.output_dir}/alevin_fry/${fastq_MD5sum}_fry_${filt_type}_quant_usa_cr-like", mode: 'copy'

    input:
        val chemistry
        val reference
        val dataset_name
        val dataset_url
        val fastq_url
        val fastq_MD5sum
        val delete_fastq
        val feature_barcode_csv_url
        val multiplexing_library_csv_url
        val filt_type
        path quant_dir

    output:
        path "dataset_description.json", emit: dataset_description

    exec:
        // from here
        // https://groups.google.com/g/nextflow/c/tp_b1p0DBE4?pli=1
        def json = JsonOutput.toJson(  [chemistry: "${chemistry}",
                                        reference: "${reference}",
                                        dataset_name: "${dataset_name}",
                                        dataset_url: "${dataset_url}",
                                        fastq_url: "${fastq_url}",
                                        fastq_MD5sum: "${fastq_MD5sum}",
                                        delete_fastq: "${delete_fastq}",
                                        feature_barcode_csv_url: "${feature_barcode_csv_url}",
                                        multiplexing_library_csv_url: "${multiplexing_library_csv_url}"
                                        ]
                                    )
        task.workDir.resolve("dataset_description.json").text = json
        // new File("$dataset_description").write(json)
}

workflow af {
    take:
        chemistry
        reference
        dataset_name
        dataset_url
        fastq_url
        fastq_MD5sum
        delete_fastq
        feature_barcode_csv_url
        multiplexing_library_csv_url
        pl_path
        t2g_path
        map_dir_path
        map_rad_path
        unmapped_file_path

    main:
        alevin_fry_gpl(
            chemistry,
            reference,
            dataset_name,
            dataset_url,
            fastq_url,
            fastq_MD5sum,
            delete_fastq,
            feature_barcode_csv_url,
            multiplexing_library_csv_url,
            "unfilt",
            pl_path,
            t2g_path,
            map_dir_path
        ) \
        | alevin_fry_collate \
        | alevin_fry_quant \
        | write_description

    emit:
        dataset_description = write_description.out.dataset_description
}


