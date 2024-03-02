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
        tuple val(chemistry),
            val(reference),
            val(dataset_name),
            val(dataset_url),
            val(fastq_url),
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
        val fastq_url, emit:fastq_url
        val fastq_MD5sum, emit: fastq_MD5sum
        val delete_fastq, emit: delete_fastq
        val feature_barcode_csv_url, emit: feature_barcode_csv_url
        val multiplexing_library_csv_url, emit: multiplexing_library_csv_url
        path index_dir_path, emit: index_dir_path
        path pl_path, emit: pl_path
        path t2g_path, emit: t2g_path
        path "${fastq_MD5sum}_fastqs", emit: fastq_dir_path

    """
        num_attempt=1
        wget ${fastq_url} -P ${fastq_MD5sum}_cfastqs
        while [ "\$(md5sum ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) | cut -d' ' -f1)" != "${fastq_MD5sum}" ]
        do
            if [[ \$num_attempt -gt 3 ]]
            then
                echo "Three attempts were made to fetch ${fastq_url.toString().lastIndexOf('/').with {fastq_url.toString().substring(it+1, fastq_url.toString().length())}}, but the MD5sum (\$(md5sum ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) | cut -d' ' -f1)) of the downloaded file \$(ls ${fastq_MD5sum}_cfastqs) didn't match the expected MD5sum ($fastq_MD5sum). Processing of this dataset will not proceed. Please check the MD5sum and internet connectivity."
                exit 1
            fi
            let "num_attempt+=1"
            rm -rf ${fastq_MD5sum}_cfastqs
            wget ${fastq_url} -P ${fastq_MD5sum}_cfastqs
        done
        mkdir -p ${fastq_MD5sum}_fastqs
        tar xf ${fastq_MD5sum}_cfastqs/\$(ls ${fastq_MD5sum}_cfastqs) --strip-components=1 -C ${fastq_MD5sum}_fastqs
        rm -rf ${fastq_MD5sum}_cfastqs
    """


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
        val fastq_url
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
        val fastq_url, emit:fastq_url
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
            reads1="\$(find -L ${fastq_dir_path} -name "*_R1_*" -type f | sort | awk '{\$1=\$1;print}' | paste -sd' ')"
            reads2="\$(find -L ${fastq_dir_path} -name "*_R2_*" -type f | sort | awk '{\$1=\$1;print}' | paste -sd' ')"
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

            reads1="\$(find -L ${fastq_dir_path} -name "*_R1_*" -type f | sort | awk '{\$1=\$1;print}' | paste -sd\" \")"
            reads2="\$(find -L ${fastq_dir_path} -name "*_R2_*" -type f | sort | awk '{\$1=\$1;print}' | paste -sd\" \")"

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


workflow salmon_map {
    take:
        samp

    main:
        download_fastq(samp) | salmon_map_rad

    emit:
        chemistry = salmon_map_rad.out.chemistry
        reference = salmon_map_rad.out.reference
        dataset_name = salmon_map_rad.out.dataset_name
        dataset_url = salmon_map_rad.out.dataset_url
        fastq_url = salmon_map_rad.out.fastq_url
        fastq_MD5sum = salmon_map_rad.out.fastq_MD5sum
        delete_fastq = salmon_map_rad.out.delete_fastq
        feature_barcode_csv_url = salmon_map_rad.out.feature_barcode_csv_url
        multiplexing_library_csv_url = salmon_map_rad.out.multiplexing_library_csv_url
        pl_path = salmon_map_rad.out.pl_path
        t2g_path = salmon_map_rad.out.t2g_path
        map_dir_path = salmon_map_rad.out.map_dir_path
        map_rad_path = salmon_map_rad.out.map_rad_path
        unmapped_file_path = salmon_map_rad.out.unmapped_file_path
}
