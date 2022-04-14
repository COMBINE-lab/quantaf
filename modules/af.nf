
/*
* This process takes the output produced by the 
* `salmon_map_rad` rule and generates a permit list 
* with alevin-fry (currently using the knee method).
*/
process alevin_fry_gpl {
    tag "fry_gpl:${MD5}:${filt_type}"
    
    conda "bioconda::alevin-fry"

    cpus 1

    input:
        val MD5
        val filt_type
        path permitlist
        path map_dir

    output:
        path "${MD5}_permitlist_${filt_type}_fw", emit: permit_dir

    script:
        filt_flag = (filt_type == "knee") ? "-k" : "-u ${permitlist}"

        """
        $params.timecmd -v -o ${MD5}_logs_permitlist_${filt_type}.time alevin-fry generate-permit-list ${filt_flag} -d fw -i ${map_dir} -o ${MD5}_permitlist_${filt_type}_fw
        """
    
    stub:
        def filt_flag = (filt_type == "knee") ? "-k" : "-u ${permitlist}"
        def cmd = """
        $params.timecmd -v -o ${MD5}_logs_permitlist_${filt_type}.time alevin-fry generate-permit-list ${filt_flag} -d fw -i ${map_dir} -o ${MD5}_permitlist_${filt_type}_fw
        """

        """
        echo "running ${cmd}"
        echo ${map_dir}
        echo ${MD5}_permitlist_${filt_type}_fw
        mkdir ${MD5}_permitlist_${filt_type}_fw
        touch ${MD5}_permitlist_${filt_type}_fw/permit_map.bin
        """
}

/*
* This process takes the output produced by the 
* `alevin_fry_gpl` rule and generates a collated 
* RAD file.
*/
process alevin_fry_collate {
    tag "fry_collate:${MD5}:${filt_type}"
    cpus params.n_threads
    conda "bioconda::alevin-fry"

    input:
        val MD5
        val filt_type
        path map_dir
        path permit_dir

    output:
        path permit_dir, emit: collate_dir

    script:
        """
        $params.timecmd -v -o ${MD5}_logs_collate_${filt_type}.time alevin-fry collate \
    -i ${permit_dir} -r ${map_dir} -t ${task.cpus} 
        """

    stub:
    def cmd = """
        $params.timecmd -v -o ${MD5}_logs_collate_${filt_type}.time alevin-fry collate \
    -i ${permit_dir} -r ${map_dir} -t ${task.cpus} 
    """

    """
    echo "executing :: ${cmd}"
    # mkdir -p ${permit_dir}
    touch ${permit_dir}/map.collated.rad
    """
}

/*
* This process takes the output produced by the 
* `alevin_fry_collate` rule and generates the resulting
* quantification.
*/
process alevin_fry_quant {
    tag "fry_quant:${MD5}:${filt_type}"
    conda "bioconda::alevin-fry"
    cpus params.n_threads

    publishDir "${params.output_dir}/alevin_fry", mode: 'rellink'
    
    input:
        val dataset_name
        val dataset_webpage
        val MD5
        val filt_type
        path t2g
        path collate_dir

    output:
        // path "${odir}/dataset_description.txt", emit: dataset_description
        // path "${odir}/**", emit: quant_dir

    script:
        odir = "${MD5}_fry_${filt_type}_quant_usa_cr-like/"
        cmd = """$params.timecmd -v -o ${MD5}_logs_quant_${filt_type}.time alevin-fry quant \
    -r cr-like --use-mtx -m ${t2g} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""

        """
            ${cmd}
            echo "${dataset_name}\n${dataset_webpage}" > ${odir}/dataset_description.txt
        """

    stub:
        odir = "${MD5}_fry_${filt_type}_quant_usa_cr-lik"
        cmd = """$params.timecmd -v -o ${MD5}_logs_quant_${filt_type}.time alevin-fry quant \
        -r cr-like --use-mtx -m ${t2g} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""

        """
            echo "executing :: ${cmd}"
            mkdir -p ${odir}
            echo "${dataset_name}\n${dataset_webpage}" > ${odir}/dataset_description.txt
            touch ${odir}/meta_info.json
        """
}

workflow af_gpl {
    take: 
        MD5
        filt_type
        permitlist
        map_dir

    main: 
        alevin_fry_gpl(
            MD5, 
            filt_type, 
            permitlist, 
            map_dir
        )

    emit: 
        alevin_fry_gpl.out.permit_dir
}

workflow af_collate {
    take: 
        MD5
        filt_type
        map_dir
        permit_dir
    
    main: 
        alevin_fry_collate(
            MD5, 
            filt_type, 
            map_dir, 
            permit_dir
        )

    emit: 
        alevin_fry_collate.out.collate_dir
}

workflow af_quant {
    take: 
        dataset_name
        dataset_webpage
        MD5
        filt_type
        t2g
        collate_dir

    main: 
        alevin_fry_quant(
            dataset_name, 
            dataset_webpage, 
            MD5, 
            filt_type, 
            t2g, 
            collate_dir
        )

    // emit: 
        // alevin_fry_quant.out.dataset_description // name
        // alevin_fry_quant.out.quant_dir // quant-dir
}

workflow af {
    take:
        dataset_name
        dataset_webpage
        MD5
        permitlist
        t2g
        map_dir
        map_rad
        unmapped_file
        filt_type

    emit:
        dataset_name

    main:
        af_gpl(MD5, filt_type, permitlist, map_dir)
        af_collate(MD5, filt_type, map_dir, af_gpl.out)
        af_quant(dataset_name, dataset_webpage, MD5, filt_type, t2g, af_collate.out)
}


