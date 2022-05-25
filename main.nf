nextflow.enable.dsl=2

include {salmon_map} from "./modules/salmon_map"
include {preprocess} from './modules/preprocess'
include {af} from './modules/af'

workflow {

      data = Channel
      .fromPath(params.input_sheets.sample)
      .splitCsv(header:true, sep: "\t", strip: true)
      .map{ row-> tuple(row.chemistry,
                        row.reference,
                        row.dataset_name,
                        row.dataset_url,
                        row.fastq_url,
                        row.fastq_MD5sum,
                        row.delete_fastq,
                        row.feature_barcode_csv_url, 
                        row.multiplexing_library_csv_url)
      }

      
      // run alevin-fry on the dataset 
      // producing both the knee filtered and 
      // unfiltered output
      preprocess()

  // merge permit list
      data = preprocess.out.chem_pl.cross(data)
            .map(it -> tuple(it[1][1], // reference
                              it[1][0], // chemistry
                              it[0][1], // pl
                              it[1][2], // dataset_name
                              it[1][3], // dataset_url
                              it[1][4], // fastq_url
                              it[1][5], // fastq_MD5sum
                              it[1][6], // delete_fastq
                              it[1][7], // feature_barcode_csv_url
                              it[1][8]  // multiplexing_library_csv_url
            ))

  // merge t2g and salmon index
      data = preprocess.out.ref_t2g_index.cross(data)
            .map(it -> tuple( it[1][1], // chemistry
                              it[0][0], // reference
                              it[1][3], // dataset_name
                              it[1][4], // dataset_url
                              it[1][5], // fastq_url
                              it[1][6], // fastq_MD5sum
                              it[1][7], // delete_fastq
                              it[1][8], // feature_barcode_csv_url
                              it[1][9], // multiplexing_library_csv_url
                              it[0][2], // index_path
                              it[1][2], // pl_path
                              it[0][1] // t2g_path
            ))

      salmon_map(data) 
      af(salmon_map.out, Channel.value("unfilt"))
}
