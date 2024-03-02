
/*
* based on the provided reference and link of permit list,
* download and maybe uncompress the files
* NOTE: require a tsv named pl_sheet.tsv in the ${bench_dir}/input_files
*/
process get_permitlist {
    tag "get_permitlist:${chemistry}"

    input:
        tuple val(chemistry), path(link)

    output:
        tuple val(chemistry), file(ofile_name), emit: chem_pl

    script:
        ofile_name = link.getName()

    if (link.getExtension() != 'txt') {
        ofile_name = link.getBaseName()
        """
        gunzip -c ${link} > ${ofile_name}
        """
    } else {
        """
        # do nothing
        """
    }
}

// turn the tuple to a map so that we can query the value by key
// permitlist = permitlist.flatMap { ["${it[0]}" : it[1]] }



/*
* based on the provided reference and link of cellranger references,
* download and uncompress on the fly, and make splici for them
* NOTE:
* 1. require a tsv named ref_sheet.tsv in the ${bench_dir}/input_files
* 2. only work for pre-built cellranger references, like humand2020A and mm10-2020A
* 3. if other reference needed, this process has to be expanded to
*    enable cellranger to make reference
*/
process get_splici {
    tag "get_permitlist:${reference}"
    // conda "bioconda::bedtools bioconda::pyroe"

    input:
        tuple val(reference), val(link)

    output:
        tuple val(reference),
        path("splici_$reference/splici_fl${params.read_len - 5}.fa"),
        path("splici_$reference/splici_fl${params.read_len - 5}_t2g_3col.tsv"),
        emit: splici

    script:
        """
        mkdir reference_$reference && wget -qO- ${link} | tar xzf - --strip-components=1 -C reference_$reference
        pyroe make-splici reference_$reference/fasta/genome.fa reference_$reference/genes/genes.gtf ${params.read_len} splici_$reference
        rm -rf reference_$reference
        """
    stub:
        """
            mkdir reference_$reference
            mkdir reference_$reference/genes && touch reference_$reference/genes/genes.gtf
            mkdir reference_$reference/fasta && touch reference_$reference/fasta/genome.fa

            mkdir splici_$reference
            touch splici_$reference/splici_fl${params.read_len-5}.fa
            touch splici_$reference/splici_fl${params.read_len-5}_t2g_3col.tsv
            rm -rf reference_$reference
        """
}


/*
* This process takes an input sample
* as defined by the columns in the sample table csv
* and runs `alevin` on it to produce a RAD file.
* Currently, sketch mode is always used.
*/
process salmon_index {
    tag "salmon_index:${reference}"

    // conda "bioconda::salmon"

    label 'multi_threads'

    input:
        tuple val(reference), path(splici), path(t2g)

    output:
        tuple val(reference), path(t2g), path("${reference}_index"), emit: ref_t2g_index

    script:

        cmd = """
            $params.timecmd -v -o ${reference}_index.time salmon index \
            -t $splici \
            -i ${reference}_index \
            -p ${task.cpus}
        """

        """
        ${cmd}
        """

    stub:
        cmd = """
            $params.timecmd -v -o ${reference}_index.time salmon index \
            -t $splici \
            -i ${reference}_index \
            -p ${task.cpus}
        """

        """
            echo "executing :: ${cmd}"
            mkdir -p ${reference}_index
            touch ${reference}_index/a_file
        """
}

/*
* This workflow take a ref_sheet.tsv and a pl_sheet.tsv file
* in the ${bench_dir}/input_files
* it prepares the permiit list and splici refernece for alevin-fry
* The output is two dictionaries, one for splici, one for permitlist.
* each item stores the name and file/folder path of the corresponding
* file/folder
*/
workflow preprocess {
    pl_sheet = Channel
            .fromPath(params.permitlist)
            .splitCsv(header:true, sep:"\t", strip: true)
            .map{ row-> tuple(row.reference,
                            row.link)
            }
    get_permitlist(pl_sheet)

    ref_sheet = Channel
                .fromPath(params.reference)
                .splitCsv(header:true, sep:"\t", strip: true)
                .map{ row-> tuple(row.reference,
                                row.link)
                }
    get_splici(ref_sheet)
    salmon_index(get_splici.out)

    chem_pl = get_permitlist.out.chem_pl
    ref_t2g_index = salmon_index.out.ref_t2g_index

    emit:
        chem_pl
        ref_t2g_index
}
