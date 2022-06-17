import groovy.io.FileType
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
* based on the provided reference and link of cellranger references, 
* download and uncompress on the fly, and make splici for them 
* NOTE: 
* 1. require a tsv named ref_sheet.tsv in the ${bench_dir}/input_files
* 2. only work for pre-built cellranger references, like humand2020A and mm10-2020A
* 3. if other reference needed, this process has to be expanded to 
*    enable cellranger to make reference 
*/
process get_splici {
    tag "get_splici:${reference}"
    // conda "bioconda::bedtools bioconda::pyroe"

    input:
        tuple val(reference), path(ref_path)
        
    output:
        tuple val(reference),
        path("splici_$reference/splici_fl${params.read_len - 5}.fa"), 
        path("splici_$reference/splici_fl${params.read_len - 5}_t2g_3col.tsv"), 
        emit: splici

    script:
        """
        pyroe make-splici ${ref_path}/fasta/genome.fa ${ref_path}/genes/genes.gtf ${params.read_len} splici_$reference
        rm -rf ${ref_path}/
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


process pre_splici {
    tag "pre_splici:${reference}"

    input:
        tuple val(reference), val(ref_data), val(fasta), val(gtf)
    
    output:
        tuple val(reference),
        path(ref_path),
        emit: processed_splici_path

    script:
        ref_path = "reference_${reference}"
        //faPath =  new File(ref_data, fasta)
        //gtfPath = new File(ref_data, gtf)
        //[ ! -f "reference_$reference$gtf" ] && echo "File not found in user specified location reference_${reference}${gtf}"
    if (ref_data.substring(0,4).equals("http")) { //If ref_data is a http link
        if (fasta.substring(0,4).equals("http") && !gtf.substring(0,4).equals("http")) { //If fasta is a https link and gtf is an absolute path
            //Download fasta, check gtf is absolute path
            gtfFile = new File(gtf)
            if (gtfFile.exits()) { //If it does exist
                unzip_file(gtf,reference)
                //if the genes directory is not made in unzip, make it
                """
                [! -d reference_$reference/genes] && mkdir reference_$reference/genes
                """
            }
        } else if (!fasta.substring(0,4).equals("http") && gtf.substring(0,4).equals("http")) { //Gtf is an https link, fasta is an absolute path
            //Download gtf, check fasta is absolute path
        } else if (fasta.substring(0,4).equals("http") && gtf.substring(0,4).equals("http")) { //Gtf and fasta are https links
            //Download both
            if(fasta.contains(".tar")) { //fasta is .tar
                if(gtf.contains(".tar")) { //gtf is .tar
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | tar xzf - --strip-components=1 -C $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | tar xzf - --strip-components=1 -C $ref_path/genes/genes.gtf
                    """
                } else if (gtf.contains(".gz")) { //gtf is .gz
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | tar xzf - --strip-components=1 -C $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | gunzip -c > $ref_path/genes/genes.gtf
                    """
                } else { //gtf is standard .gtf file
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | tar xzf - --strip-components=1 -C $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qP $gtf $ref_path/genes/genes.gtf
                    """
                }
            } else if (fasta.contains(".gz")) { //fasta is .gz
                if (gtf.contains(".gz")) { //gtf is .gz
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | gunzip -c > $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | gunzip -c > $ref_path/genes/genes.gtf
                    """
                } else if (gtf.contains(".tar")) { //gtf is tar
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | gunzip -c > $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | tar xzf - --strip-components=1 -C $ref_path/genes/genes.gtf
                    """
                } else { //gtf is standard
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    wget -qO- $fasta | gunzip -c > $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | gunzip -c > $ref_path/genes/genes.gtf
                    """
                }
            } else { //fasta is .fa
                if (gtf.contains(".gz")) { //gtf is .gz
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    mv $fasta_p $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    gunzip -c $gtf_p > $ref_path/genes/genes.gtf
                    """
                } else if (gtf.contains(".tar")) { //gtf is .tar
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    mv $fasta_p $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qO- $gtf | tar xzf - --strip-components=1 -C $ref_path/genes/genes.gtf
                    """
                } else { //gtf is .gtf
                    """
                    mkdir $ref_path
                    mkdir $ref_path/fasta/
                    mv $fasta_p $ref_path/fasta/genome.fa
                    mkdir $ref_path/genes/
                    wget -qP $gtf $ref_path/genes/genes.gtf
                    """
                }
            }
        } else if (!fasta.substring(0,4).equals("http") && !gtf.substring(0,4).equals("http")) {
            //check if its an absolute path, otherwise download ref_data and check in specified locations
            //ref_path = file("reference_$reference/")    
            """
            mkdir reference_$reference
            wget -qO- ${ref_data} | tar xzf - --strip-components=1 -C reference_$reference
            """
        }
    } else if(ref_data.substring(ref_data.length()-7, ref_data.length()).equals(".tar.gz")) {
        """
        mkdir reference_$reference
        tar xvf ${ref_data} --strip-components=1 -C reference_$reference
        """
    } else if(ref_data.substring(ref_data.length()-4, ref_data.length()).equals(".tar")) {
        """
        mkdir reference_$reference
        tar xvf ${ref_data} --strip-components=1 -C reference_$reference
        """
    } else if(ref_data.substring(ref_data.length()-3, ref_data.length()).equals(".gz")) {
        """
        mkdir reference_$reference
        gunzip -c ${ref_data} > reference_$reference
        """
    } else {
        """
        echo Invalid ref_data provided for reference $reference
        """
    }
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
            .fromPath(params.input_sheets.permitlist)
            .splitCsv(header:true, sep:"\t", strip: true)
            .map{ row-> tuple(row.chemistry,
                            row.link)
            }
    get_permitlist(pl_sheet)
    
    ref_sheet = Channel
                .fromPath(params.input_sheets.reference)
                .splitCsv(header:true, sep:"\t", strip: true)
                .map{ row-> tuple(row.reference,
                                row.ref_data, row.fasta, row.gtf)
                }
    pre_splici(ref_sheet)
    get_splici(pre_splici.out)
    salmon_index(get_splici.out)
    
    chem_pl = get_permitlist.out.chem_pl
    ref_t2g_index = salmon_index.out.ref_t2g_index
    emit:
        chem_pl
        ref_t2g_index
}