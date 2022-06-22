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
        tuple val(reference), path(ref_path), path(fasta), path(gtf)
        
    output:
        tuple val(reference),
        path("splici_$reference/splici_fl${params.read_len - 5}.fa"), 
        path("splici_$reference/splici_fl${params.read_len - 5}_t2g_3col.tsv"), 
        emit: splici

    script:
        """
        pyroe make-splici $fasta $gtf ${params.read_len} splici_$reference
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

process process_ref_data {
    tag "ref:${reference}"

    input: 
        tuple val(reference), val(ref_data), val(fasta), val(gtf)
    output:
        tuple val(reference), path(ref_path), val(fasta), val(gtf),
        emit: processed_ref_data        

    script:
        ref_path = "reference_${reference}/"
    if (!fasta.substring(0,4).equals("http") && !gtf.substring(0,4).equals("http")) {   
        if (ref_data.substring(0,4).equals("http")) { //check if ref_data is https or local
            """
            mkdir $ref_path
            wget -qO- $ref_data | tar xzf - --strip-components=1 -C $ref_path
            find > ~/out2.txt
            pwd >> ~/out2.txt
            """
        } else { //if its local see if its a tarball or if its just fasta and gtf
            if(ref_data.contains("gz") || ref_data.contains("tar")) {
                """
                mkdir $ref_path
                tar xzf - --strip-components=1 -C $ref_path
                """
            } else {
                """
                mkdir $ref_path
                """
                //cp -rf $ref_data $ref_path/
            }    
        }
    } else {
        """
        mkdir $ref_path
        """
    }  
}

process standardize_files {
    tag "standardize_files:$reference"

    input: 
        tuple val(reference), path(ref_path), val(fasta), val(gtf)
    output:
        tuple val(reference), path(ref_path), val(fasta), val(gtf),
        emit: standardized_files_data
    script:
        File abs_fasta = new File("$fasta") //absolute path to given fasta
        File rel_fasta = new File("${projectDir}/${fasta}") //path relative to project dir given for fasta
        File ref_fasta = new File("${ref_path}${fasta}")  //path relative to reference directory given for fasta
        File abs_gtf = new File(gtf) //absolute path to given gtf
        File rel_gtf= new File("${projectDir}/${gtf}") //path relative to project dir given for gtf
        File ref_gtf = new File("${ref_path}${gtf}")  //path relative to reference directory given for gtf    
        if (rel_fasta.exists()) {
            fasta = "${projectDir}/${fasta}"
        } else if (ref_fasta.exists()) {
            fasta = "${ref_path}${fasta}"
        } else {
            //do nothing for now
        }
        if (rel_gtf.exists()) {
            gtf = "${projectDir}/${gtf}"
        } else if (ref_gtf.exists()) {
            gtf = "${ref_path}/${gtf}"
        }
        """
        pwd > ~/out.txt
        ls -al >> ~/out.txt
        find >> ~/out.txt
        """
        
}

process process_fasta {
    tag "process_fasta:$fasta"

    input:
        tuple val(reference), path(ref_path), val(fasta), val(gtf)
    output:
        tuple val(reference), path(ref_path), path(fasta_p), val(gtf),
        emit: processed_fasta_path    
    script:
        File abs_fasta = new File(fasta)
        fasta_p = "$ref_path/fasta/genome.fa"
        if (fasta.substring(0,4).equals("http") && fasta.contains(".tar")) { //fasta is remote tarball
            //download and untar
            """
            mkdir $ref_path/fasta/
            wget -qO- $fasta | tar xzf - --strip-components=1 -C $ref_path/fasta/genome.fa
            """
        } else if (fasta.substring(0,4).equals("http") && fasta.contains(".gz")) { //fasta is remote gzip
            //download and unzip
            """
            mkdir $ref_path/fasta/
            wget -qO- $fasta | gunzip -c > $ref_path/fasta/genome.fa
            """
        }  else if (fasta.substring(0,4).equals("http")) { //fasta is remote
            //download and move
            """
            mkdir $ref_path/fasta/
            wget -qP $fasta $ref_path/fasta/genome.fa
            """
        } else if (abs_fasta.exists()) { //if the path given is an absolute path to a file
            if (fasta.contains(".tar")) { //tarball 
                //unta
                """
                mkdir $ref_path/fasta/
                tar -xzf $fasta --strip-components=1 -C $ref_path/fasta/genome.fa
                """
            } else if (fasta.contains(".gz")) { //gzip
                //gunzip
                """
                mkdir $ref_path/fasta/
                gunzip $fasta -c > $ref_path/fasta/genome.fa
                """
            } else { //.fa
                //move file to location exisiting
                """
                mkdir $ref_path/fasta/
                cp $fasta $ref_path/fasta/genome.fa
                """
            }
        }
}

process process_gtf {
    tag "process_gtf:$gtf"

    input:
        tuple val(reference), path(ref_path), path(fasta), val(gtf)
    output:
        tuple val(reference), path(ref_path), path(fasta_p), path(gtf_p),
        emit: processed_gtf_path    
    script:
        File abs_gtf = new File(gtf)
        fasta_p = fasta
        gtf_p = "$ref_path/genes/genes.gtf"
        if (gtf.substring(0,4).equals("http") && gtf.contains(".tar")) { //gtf is remote tarball
            //download and untar
            """
            mkdir $ref_path/genes/
            wget -qO- $gtf | tar xzf - --strip-components=1 -C $ref_path/genes/genes.gtf
            """
        } else if (gtf.substring(0,4).equals("http") && gtf.contains(".gz")) { //gtf is remote gzip
            //download and unzip
            """
            mkdir $ref_path/genes/
            wget -qO- $gtf | gunzip -c > $ref_path/genes/genes.gtf
            """
        }  else if (gtf.substring(0,4).equals("http")) { //fasta is remote
            //download and move
            """
            mkdir $ref_path/genes/
            wget -qP $gtf $ref_path/genes/genes.gtf
            """
        } else if (abs_fasta.exists()) { //if the path given is an absolute path to a file
            if (gtf.contains(".tar")) { //tarball 
                //untar
                """
                mkdir $ref_path/genes/
                tar -xzf $gtf --strip-components=1 -C $ref_path/genes/genes.gtf
                """
            } else if (fasta.contains(".gz")) { //gzip
                //gunzip
                """
                mkdir $ref_path/genes/
                gunzip $gtf -c > $ref_path/genes/genes.gtf
                """
            } else { //.fa
                //move file to location exisiting
                """
                mkdir $ref_path/genes/
                cp $gtf $ref_path/genes/genes.gtf
                """
            }
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
    process_ref_data(ref_sheet)
    standardize_files(process_ref_data.out)

    process_fasta(standardize_files.out)
    process_gtf(process_fasta.out)

    get_splici(process_gtf.out)         
    salmon_index(get_splici.out)
    
    chem_pl = get_permitlist.out.chem_pl
    ref_t2g_index = salmon_index.out.ref_t2g_index
    emit:
        chem_pl
        ref_t2g_index
}