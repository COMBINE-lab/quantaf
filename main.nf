nextflow.enable.dsl=2

// include { salmon_map; af_knee; af_unfilt } from './modules/alevin'
include {salmon_map} from "./modules/salmon_map"
include {preprocess} from './modules/preprocess'
include {af} from './modules/af'

workflow {

    data = Channel
    .fromPath(params.input_sheets.sample)
    .splitCsv(header:true, sep: "\t")
    .map{ row-> tuple(row.chemistry,
                      row.reference,
                      row.dataset_name,
                      row.link,
                      row.data_url,
                      row.MD5,
                      row.delete_fastq,
                      row.feature_barcode, 
                      row.library_csv)
    }

  
  // run alevin-fry on the dataset 
  // producing both the knee filtered and 
  // unfiltered output
  preprocess()

  // merge permit list
  data = preprocess.out.chem_pl.cross(data)
        .map(it -> tuple(it[1][1], // ref
                          it[1][0], // chem
                          it[0][1], // pl
                          it[1][2], // dataset_name
                          it[1][3], // dataset_webpage
                          it[1][4], // fastq_link
                          it[1][5], // MD5
                          it[1][6], // delete_fastq
                          it[1][7], // feature_barcode
                          it[1][8]  // multiplexing
        ))
  
  // merge t2g and salmon index
  data = preprocess.out.ref_t2g_index.cross(data)
        .map(it -> tuple(it[1][3], // dataset_name
                          it[1][5], // fastq_link
                          it[1][6], // MD5
                          it[1][7], // delete_fastq
                          it[1][1], // chem
                          it[0][2], // index
                          it[1][2], // pl
                          it[0][1], // t2g
                          it[0][0], // ref
                          it[1][4], // dataset_webpage
                          it[1][8], // feature_barcode
                          it[1][9] // multiplexing
        ))

  salmon_map(data)
  af(salmon_map.out, Channel.value("unfilt"))
}
