workflow salmon_map {
    take: 
        samp 

    main: 
        salmon_map_rad(samp)

    emit:
        salmon_map_rad.out.dataset_name
        salmon_map_rad.out.dataset_webpage
        salmon_map_rad.out.MD5
        salmon_map_rad.out.permitlist
        salmon_map_rad.out.t2g
        salmon_map_rad.out.map_dir
        salmon_map_rad.out.map_rad
        salmon_map_rad.out.unmapped_file
}

/*
* This process takes an input sample
* as defined by the columns in the sample table csv
* and runs `alevin` on it to produce a RAD file.
* Currently, sketch mode is always used.
*/
process salmon_map_rad {
    tag "salmon_map:${MD5}"

    label 'multi_threads'

    conda "bioconda::salmon"

    input:
        tuple val(dataset_name), 
            val(fastq_link), 
            val(MD5), 
            val(delete_fastq), 
            val(chem), 
            path(index), 
            path(pl), 
            path(t2g), 
            val(ref), 
            val(dataset_webpage), 
            val(feature_barcode), 
            val(multiplexing)
    output:
        val dataset_name, emit: dataset_name
        val dataset_webpage, emit: dataset_webpage
        val MD5, emit:MD5
        path pl, emit: permitlist
        path t2g, emit: t2g
        path "${MD5}_alevin_map", emit: map_dir
        path "${MD5}_alevin_map/map.rad", emit: map_rad
        path "${MD5}_alevin_map/unmapped_bc_count.bin", emit: unmapped_file

    script:
        chemistry = salmon_chem_flag(chem)

        """
            wget ${fastq_link} -P ${MD5}_cfastqs
            while [ "\$(md5sum ${MD5}_cfastqs/\$(ls ${MD5}_cfastqs) | cut -d' ' -f1)" != "${MD5}" ]
            do
                rm -rf ${MD5}_cfastqs
                wget ${fastq_link} -P ${MD5}_cfastqs
            done
            mkdir -p ${MD5}_fastqs
            tar xf ${MD5}_cfastqs/\$(ls ${MD5}_cfastqs) --strip-components=1 -C ${MD5}_fastqs
            rm -rf ${MD5}_cfastqs

            reads1="\$(find ${MD5}_fastqs -name "*_R1_*" -type f | xargs | sort | awk '{print " "\$0}')"
            reads2="\$(find ${MD5}_fastqs -name "*_R2_*" -type f | xargs | sort | awk '{print " "\$0}')"

            /usr/bin/time -v -o ${MD5}_log_map_sketch.time \
            salmon alevin -i ${index} -l ISR \
            -1 \$reads1 \
            -2 \$reads2 \
            -p ${task.cpus} \
            --${chemistry} \
            --sketch -o \
            ${MD5}_alevin_map 

            if [ ${delete_fastq != 0} ]; then rm -rf ${MD5}_fastqs; fi
        """

    stub:
        chemistry = salmon_chem_flag(chem)

        """
            echo ${chemistry}
            mkdir -p ${MD5}_fastqs
            touch ${MD5}_fastqs/read_S1_L001_R1_001.fastq.gz
            touch ${MD5}_fastqs/read_S1_L001_R2_001.fastq.gz
            echo "\$(ls ${MD5}_fastqs)"

            reads1="\$(ls ${MD5}_fastqs | sort | awk -v p=${MD5}_fastqs '{print p"/"\$0}' | grep "_R1_")"
            reads2="\$(ls ${MD5}_fastqs | sort | awk -v p=${MD5}_fastqs '{print p"/"\$0}' | grep "_R2_")"

            echo \$reads1
            echo \$reads2

            mkdir -p ${MD5}_alevin_map
            touch ${MD5}_alevin_map/map.rad
            touch ${MD5}_alevin_map/unmapped_bc_count.bin

            echo "\$(ls ${MD5}_alevin_map)"

            if [ ${delete_fastq != 0} ]; then rm -rf ${MD5}_fastqs; fi
        """

}


/*
* extract the flag that should be passed to alevin
* for the chemistry based on the chemistry type of a
* sample.  Currently, everything is Chromium v2 or v3
* and so that is all this currently supports
*/
def salmon_chem_flag(chem) {
    switch(chem) {
        case "v2":
        return "chromium"
        case "v3":
        return "chromiumV3"
    }
}

