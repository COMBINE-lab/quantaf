
/*
* This process takes the output produced by the 
* `salmon_map_rad` rule and generates a permit list 
* with alevin-fry (currently using the knee method).
*/
process alevin_fry_gpl {
    tag "fry_gpl:${MD5}:${filt_type}"
    conda "bioconda::alevin-fry"
    label 'single_threads'

    input:
        val dataset_name
        val dataset_webpage
        val MD5
        path permitlist
        path t2g
        path map_dir
        val filt_type

    output:
        val dataset_name, emit: dataset_name
        val dataset_webpage, emit: dataset_webpage
        val MD5, emit: MD5
        path permitlist, emit: permitlist
        path t2g, emit: t2g
        path map_dir, emit: map_dir
        val filt_type, emit: filt_type
        path "${MD5}_permitlist_${filt_type}_fw", emit: permit_dir

    script:
        filt_flag = (filt_type == "knee") ? "-k" : "-u ${permitlist}"

        cmd = """
        $params.timecmd -v -o ${MD5}_logs_permitlist_${filt_type}.time alevin-fry generate-permit-list ${filt_flag} -d fw -i ${map_dir} -o ${MD5}_permitlist_${filt_type}_fw
        """

        """
        ${cmd}
        """
    
    stub:
        filt_flag = (filt_type == "knee") ? "-k" : "-u ${permitlist}"
        cmd = """
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
    label 'multi_threads'
    conda "bioconda::alevin-fry"

    input:
        val dataset_name
        val dataset_webpage
        val MD5
        path permitlist
        path t2g
        path map_dir
        val filt_type
        path permit_dir

    output:
        val dataset_name, emit: dataset_name
        val dataset_webpage, emit: dataset_webpage
        val MD5, emit: MD5
        path permitlist, emit: permitlist
        path t2g, emit: t2g
        path map_dir, emit: map_dir
        val filt_type, emit: filt_type
        path permit_dir, emit: collate_dir

    script:
        cmd = """
        $params.timecmd -v -o ${MD5}_logs_collate_${filt_type}.time alevin-fry collate \
    -i ${permit_dir} -r ${map_dir} -t ${task.cpus} 
        """

        """
        ${cmd}
        """

    stub:
    cmd = """
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
    label 'multi_threads'

    publishDir "${params.output_dir}/alevin_fry", mode: 'copy'
    
    input:
        val dataset_name
        val dataset_webpage
        val MD5
        path permitlist
        path t2g
        path map_dir
        val filt_type
        path collate_dir

    output:
        path "${odir}/dataset_description.txt", emit: dataset_description
        path "${odir}/**", emit: quant_dir

    script:
        odir = "${MD5}_fry_${filt_type}_quant_usa_cr-like"
        cmd = """$params.timecmd -v -o ${MD5}_logs_quant_${filt_type}.time alevin-fry quant \
    -r cr-like --use-mtx -m ${t2g} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""

        """
            ${cmd}
            echo ${dataset_name} > ${odir}/dataset_description.txt
            echo "${dataset_webpage}" >> ${odir}/dataset_description.txt
            tar cf ${odir}.tar ${odir}
        """

    stub:
        odir = "${MD5}_fry_${filt_type}_quant_usa_cr-like"
        cmd = """$params.timecmd -v -o ${MD5}_logs_quant_${filt_type}.time alevin-fry quant \
    -r cr-like --use-mtx -m ${t2g} -i ${collate_dir} -o ${odir} -t ${task.cpus}"""

        """
            echo "executing :: ${cmd}"
            mkdir -p ${odir}
            touch ${odir}/dataset_description.txt
            echo ${dataset_name} > ${odir}/dataset_description.txt
            echo "${dataset_webpage}" >> ${odir}/dataset_description.txt
            touch ${odir}/meta_info.json
            tar cf ${odir}.tar ${odir}

        """
}

workflow af_gpl {
    take: 
        dataset_name
        dataset_webpage
        MD5
        permitlist
        t2g
        map_dir
        filt_type


    main: 
        alevin_fry_gpl(
            dataset_name,
            dataset_webpage,
            MD5,
            permitlist,
            t2g,
            map_dir,
            filt_type
        )

    emit: 
        alevin_fry_gpl.out.dataset_name
        alevin_fry_gpl.out.dataset_webpage
        alevin_fry_gpl.out.MD5
        alevin_fry_gpl.out.permitlist
        alevin_fry_gpl.out.t2g
        alevin_fry_gpl.out.map_dir
        alevin_fry_gpl.out.filt_type
        alevin_fry_gpl.out.permit_dir
}

workflow af_collate {
    take: 
        dataset_name
        dataset_webpage
        MD5
        permitlist
        t2g
        map_dir
        filt_type
        permit_dir
    
    main: 
        alevin_fry_collate(
            dataset_name,
            dataset_webpage,
            MD5,
            permitlist,
            t2g,
            map_dir,
            filt_type, 
            permit_dir
        )

    emit: 
        alevin_fry_collate.out.dataset_name
        alevin_fry_collate.out.dataset_webpage
        alevin_fry_collate.out.MD5
        alevin_fry_collate.out.permitlist
        alevin_fry_collate.out.t2g
        alevin_fry_collate.out.map_dir
        alevin_fry_collate.out.filt_type
        alevin_fry_collate.out.collate_dir
}

workflow af_quant {
    take: 
        dataset_name
        dataset_webpage
        MD5
        permitlist
        t2g
        map_dir
        filt_type
        collate_dir

    main: 
        alevin_fry_quant(
            dataset_name,
            dataset_webpage,
            MD5,
            permitlist,
            t2g,
            map_dir,
            filt_type, 
            collate_dir
        )

    emit: 
        // alevin_fry_quant.out.dataset_description // name
        alevin_fry_quant.out.quant_dir // quant-dir
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
        af_gpl(
            dataset_name,
            dataset_webpage,
            MD5,
            permitlist,
            t2g,
            map_dir,
            filt_type
        ) \
        | af_collate \
        | af_quant
}


