workflow salmon_map {
    take: 
        samp 

    main: 
        standardize_salmon_files(samp) | download_fastq | salmon_map_rad

    emit:
        salmon_map_rad.out.chemistry
        salmon_map_rad.out.reference
        salmon_map_rad.out.dataset_name
        salmon_map_rad.out.dataset_url
        salmon_map_rad.out.fastq_file
        salmon_map_rad.out.fastq_MD5sum
        salmon_map_rad.out.delete_fastq
        salmon_map_rad.out.feature_barcode_csv_url
        salmon_map_rad.out.multiplexing_library_csv_url
        salmon_map_rad.out.pl_path
        salmon_map_rad.out.t2g_path
        salmon_map_rad.out.map_dir_path
        salmon_map_rad.out.map_rad_path
        salmon_map_rad.out.unmapped_file_path

}

process standardize_salmon_files {
    tag "standardize_salmon_files:${dataset_name}"
    errorStrategy 'terminate'
    
    input:
        tuple val(chemistry), 
            val(reference), 
            val(dataset_name), 
            val(dataset_url), 
            val(fastq_file), 
            val(fastq_MD5sum), 
            val(delete_fastq), 
            val(feature_barcode_csv_url), 
            val(multiplexing_library_csv_url),
            path(index_dir_path), 
            path(pl_path), 
            path(t2g_path)
    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_file, emit:fastq_file
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        path index_dir_path, emit: index_dir_path
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path

    script:
        if(fastq_file.length () > 3 && !fastq_file.substring(0,4).equals("http")) { //fastq is not a url, and is a file
            if(!fastq_file.substring(0,1).equals("/")) {
                fastq_file = "/" + "$fastq_file"
            }
            rel_fastq = file("${projectDir}${fastq_file}") //fastq is relative to project directory
            abs_fastq = file("$fastq_file") //fastq is absolute path
            print(rel_fastq)
            if (rel_fastq.exists()) { //fastq exists relative to project directory
                fastq_file = rel_fastq.toRealPath()
                print("HERe")
            } else if (!abs_fastq.exists()) {
                error "Could not find referenced fastq file."
            }
        }
        //todo csvs for multiplex and feature barcodes?

        """

        """
    
}

/*
* This process takes an input sample
* as defined by the columns in the sample table csv
* and runs `alevin` on it to produce a RAD file.
* Currently, sketch mode is always used.
*/

process download_fastq {
    tag "download_fastq:${fastq_MD5sum}"
    label 'single_threads'

    input:
        val chemistry 
        val reference
        val dataset_name 
        val dataset_url 
        val fastq_file 
        val fastq_MD5sum 
        val delete_fastq 
        val feature_barcode_csv_url 
        val multiplexing_library_csv_url
        path index_dir_path
        path pl_path
        path t2g_path

    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_file, emit:fastq_file
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        path index_dir_path, emit: index_dir_path
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path
        path "${fastq_MD5sum}_fastqs", emit: fastq_dir_path
    
    script:
        fastq_file = fastq_file.toString()
        abs_fastq = file("$fastq_file")
        if(fastq_file.length() > 3 && fastq_file.substring(0,4).equals("http")) { //fastq is a URL
            """
                num_attempt=1
                wget ${fastq_file} -P ${fastq_MD5sum}_cfastqs
                while [ "\$(md5sum ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) | cut -d' ' -f1)" != "${fastq_MD5sum}" ]
                do
                    if [[ \$num_attempt -gt 3 ]]
                    then
                        echo "Three attempts were made to fetch ${fastq_file.toString().lastIndexOf('/').with {fastq_file.toString().substring(it+1, fastq_file.toString().length())}}, but the MD5sum (\$(md5sum ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) | cut -d' ' -f1)) of the downloaded file \$(ls ${fastq_MD5sum}_cfastqs) didn't match the expected MD5sum ($fastq_MD5sum). Processing of this dataset will not proceed. Please check the MD5sum and internet connectivity."
                        exit 1
                    fi
                    let "num_attempt+=1"
                    rm -rf ${fastq_MD5sum}_cfastqs
                    wget ${fastq_file} -P ${fastq_MD5sum}_cfastqs
                done
                mkdir -p ${fastq_MD5sum}_fastqs
                tar xf ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) --strip-components=1 -C ${fastq_MD5sum}_fastqs
                rm -rf ${fastq_MD5sum}_cfastqs
            """
        } else if(abs_fastq.exists()) { //check that the fastq file exists locally
            if (fastq_MD5sum.length() != 32) {
                fastq_MD5sum = fastq_file.md5()
            }
            if(fastq_file.length () > 3 && fastq_file.contains(".tar") && fastq_file.contains(".gz")) {
                """
                mkdir ${fastq_MD5sum}_fastqs
                tar -xzf ${fastq_file} --strip-components=1 -C ${fastq_MD5sum}_fastqs
                """
            } else if(fastq_file.length () > 3 && fastq_file.contains(".tar")) {
                """
                mkdir ${fastq_MD5sum}_fastqs
                tar -xf ${fastq_file} --strip-components=1 -C ${fastq_MD5sum}_fastqs
                """
            } else if (fastq_file.length () > 3 && fastq_file.contains(".gz")) {
                """
                mkdir ${fastq_MD5sum}_fastqs
                gunzip ${fastq_file} -c > ${fastq_MD5sum}_fastqs
                """
            } else {
                """
                mkdir ${fastq_MD5sum}_fastqs
                cp -a ${fastq_file} ${fastq_MD5sum}_fastqs
                """
            }
        }
}

process salmon_map_rad {
    tag "salmon_map:${fastq_MD5sum}"
    label 'multi_threads'
    // conda "bioconda::salmon"
    
    input:
        val chemistry
        val reference
        val dataset_name
        val dataset_url
        val fastq_file
        val fastq_MD5sum
        val delete_fastq
        val feature_barcode_csv_url
        val multiplexing_library_csv_url
        path index_dir_path
        path pl_path
        path t2g_path
        path fastq_dir_path

    output:
        val chemistry, emit: chemistry
        val reference, emit: reference
        val dataset_name, emit: dataset_name
        val dataset_url, emit: dataset_url
        val fastq_file, emit:fastq_file
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path
        path "${fastq_MD5sum}_alevin_map", emit: map_dir_path
        path "${fastq_MD5sum}_alevin_map/map.rad", emit: map_rad_path
        path "${fastq_MD5sum}_alevin_map/unmapped_bc_count.bin", emit: unmapped_file_path

    script:
        chemistry_salmon = salmon_chem_flag(chemistry)

        """
            reads1="\$(find -L ${fastq_dir_path} -name "*_R1_*" -type f | xargs | sort | awk '{print \$0}')"
            reads2="\$(find -L ${fastq_dir_path} -name "*_R2_*" -type f | xargs | sort | awk '{print \$0}')"
            ${params.timecmd} -v -o ${fastq_MD5sum}_log_map_sketch.time \
            salmon alevin -i ${index_dir_path} -l ISR \
            -1 \${reads1} \
            -2 \${reads2} \
            -p ${task.cpus} \
            --${chemistry_salmon} \
            --sketch -o \
            ${fastq_MD5sum}_alevin_map 

            if [ ${delete_fastq != 0} ]; then rm -rf ${fastq_MD5sum}_fastqs; fi
        """

    stub:
        chemistry_salmon = salmon_chem_flag(chemistry)

        """
            echo ${chemistry_salmon}
            mkdir -p ${fastq_MD5sum}_fastqs
            touch ${fastq_MD5sum}_fastqs/read_S1_L001_R1_001.fastq.gz
            touch ${fastq_MD5sum}_fastqs/read_S1_L001_R2_001.fastq.gz
            echo "\$(ls ${fastq_MD5sum}_fastqs)"

            reads1="\$(ls ${fastq_MD5sum}_fastqs | sort | awk -v p=${fastq_MD5sum}_fastqs '{print p"/"\$0}' | grep "_R1_")"
            reads2="\$(ls ${fastq_MD5sum}_fastqs | sort | awk -v p=${fastq_MD5sum}_fastqs '{print p"/"\$0}' | grep "_R2_")"

            echo \$reads1
            echo \$reads2

            mkdir -p ${fastq_MD5sum}_alevin_map
            touch ${fastq_MD5sum}_alevin_map/map.rad
            touch ${fastq_MD5sum}_alevin_map/unmapped_bc_count.bin

            echo "\$(ls ${fastq_MD5sum}_alevin_map)"

            if [ ${delete_fastq != 0} ]; then rm -rf ${fastq_MD5sum}_fastqs; fi
        """

}


/*
* extract the flag that should be passed to alevin
* for the chemistry based on the chemistry type of a
* sample.  Currently, everything is Chromium v2 or v3
* and so that is all this currently supports
*/
def salmon_chem_flag(chemistry) {
    switch(chemistry) {
        case "v2":
        return "chromium"
        case "v3":
        return "chromiumV3"
    }
}

