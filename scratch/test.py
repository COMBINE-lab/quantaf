import glob
import json
import os
import shutil
import subprocess
from pathlib import Path

from latch.types.directory import LatchDir
from latch_cli.extras.nextflow.file_persistence import (
    _extract_paths,
    download_files,
    upload_files,
)

# channel_vals = [[
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "500 Human PBMCs, 3' LT v3.1, Chromium Controller"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-controller-3-1-low-6-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_Controller/500_PBMC_3p_LT_Chromium_Controller_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "5f080c6082f11ea9fc6448482e6fb590"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "500 Human PBMCs, 3' LT v3.1, Chromium X"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "5b36a7bfda36a7093adc8e30c3fa92c8"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "1k PBMCs from a Healthy Donor (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "265ebe8f77ad90db350984d9c7a59e52"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "10k PBMCs from a Healthy Donor (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "e0021592e209642d71f5dc420cf4c5c0"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "10k Human PBMCs, 3' v3.1, Chromium X"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-ht-v3-1-chromium-x-3-1-high"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "43ea77ed6f860597c568e9ff819f0504"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "10k Human PBMCs, 3' v3.1, Chromium Controller"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "05cccabbad1b94c83f67c67b15b07e7e"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "10k Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor, Single Indexed"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-single-indexed-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_PBMC_10K/SC3_v3_NextGem_SI_PBMC_10K_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "38d2d253f8537d3c39aaa5832da39ad3"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "10k Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor, Dual Indexed"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-dual-indexed-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_PBMC_10K/SC3_v3_NextGem_DI_PBMC_10K_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "bd61d65af63508368c57f137090f32b9"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "20k Human PBMCs, 3' HT v3.1, Chromium X"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.0/20k_PBMC_3p_HT_nextgem_Chromium_X/20k_PBMC_3p_HT_nextgem_Chromium_X_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "650201ca6bd1af3f984a629da318d23b"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "PBMCs from EDTA-Treated Blood Collection Tubes Isolated via"
#                         " SepMate-Ficoll Gradient (3' v3.1 Chemistry)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbmcs-3p_edta_sepmate-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_EDTA_SepMate/3p_EDTA_SepMate_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "d07c96c9f224d75967759db6be54d635"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "PBMCs from Heparin-Treated Blood Collection Tubes Isolated via"
#                         " SepMate-Ficoll Gradient (3' v3.1 Chemistry)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbmcs-3p_heparin_sepmate-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_Heparin_SepMate/3p_Heparin_SepMate_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "4f9d8c6052713a373a6b05180dc7a1f3"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "PBMCs from ACD-A Treated Blood Collection Tubes Isolated via"
#                         " SepMate-Ficoll Gradient (3' v3.1 Chemistry)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbmcs-3p_acda_sepmate-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_ACDA_SepMate/3p_ACDA_SepMate_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "c84db3199c36884044ea4f39e8ade802"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "PBMCs from Citrate-Treated Blood Collection Tubes Isolated via"
#                         " SepMate-Ficoll Gradient (3' v3.1 Chemistry)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_sepmate-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_Citrate_SepMate/3p_Citrate_SepMate_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "367649452979be8828314810952ef551"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "PBMCs from Citrate-Treated Cell Preparation Tubes (3' v3.1"
#                         " Chemistry)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbmcs-3p_citrate_cpt-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_Citrate_CPT/3p_Citrate_CPT_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "196d59f2004579289728269016630544"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": "PBMCs from a Healthy Donor: Whole Transcriptome Analysis"
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/pbm-cs-from-a-healthy-donor-whole-transcriptome-analysis-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "5e97b62f74cca767c3c6088b54fee607"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Whole Blood RBC Lysis for PBMCs and Neutrophils,"
#                         " Granulocytes, 3'"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/WB_Lysis_Granulocytes_3p_Introns_8kCells/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "e4f2823c531e3eb289161c979c75f051"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor - Manual (channel 5)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-5-3-1-standard-3-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.1.0/manual_5k_pbmc_NGSC3_ch5/manual_5k_pbmc_NGSC3_ch5_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "2f85c73fe9c3e8a2d76064640cbe7c89"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor - Manual (channel 1)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-1-3-1-standard-3-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.1.0/manual_5k_pbmc_NGSC3_ch1/manual_5k_pbmc_NGSC3_ch1_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "ca2dd76a6a296180f23c0a782d02691f"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor - Chromium Connect (channel 5)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-5-3-1-standard-3-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.1.0/connect_5k_pbmc_NGSC3_ch5/connect_5k_pbmc_NGSC3_ch5_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "97cb61b31732ba2e6d1974b24bf531bb"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Peripheral blood mononuclear cells (PBMCs) from a healthy"
#                         " donor - Chromium Connect (channel 1)"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-chromium-connect-channel-1-3-1-standard-3-1-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.1.0/connect_5k_pbmc_NGSC3_ch1/connect_5k_pbmc_NGSC3_ch1_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "8a3973b93a26d8d4a68bddcfc1101ff8"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome"
#                         " Analysis"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/hodgkins-lymphoma-dissociated-tumor-whole-transcriptome-analysis-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/Parent_NGSC3_DI_HodgkinsLymphoma/Parent_NGSC3_DI_HodgkinsLymphoma_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "c083236af78d5556a88b3878baa40afe"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "200 Sorted Cells from Human Glioblastoma Multiforme, 3' LT"
#                         " v3.1"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/200-sorted-cells-from-human-glioblastoma-multiforme-3-lt-v-3-1-3-1-low-6-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p_LT/Brain_Tumor_3p_LT_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "c1d5f7a04b3615a5bc6cb3d3ed2e0f0a"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "750 Sorted Cells from Human Invasive Ductal Carcinoma, 3' LT"
#                         " v3.1"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/750-sorted-cells-from-human-invasive-ductal-carcinoma-3-lt-v-3-1-3-1-low-6-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Breast_Cancer_3p_LT/Breast_Cancer_3p_LT_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "b105c24862021e822982f8a1136888e2"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "2k Sorted Cells from Human Glioblastoma Multiforme, 3' v3.1"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "6bf6e6ff6b8b96c1967ef737adb5e807"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "7.5k Sorted Cells from Human Invasive Ductal Carcinoma, 3'"
#                         " v3.1"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/7-5-k-sorted-cells-from-human-invasive-ductal-carcinoma-3-v-3-1-3-1-standard-6-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.0.0/Breast_Cancer_3p/Breast_Cancer_3p_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "70a14452298a4ebf9fef72c03bafb764"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "Human Glioblastoma Multiforme: 3' v3 Whole Transcriptome"
#                         " Analysis"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "b32e92d47da181c99de5c95107d32acb"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "1k Brain Cells from an E18 Mouse (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "b74fd4a8fb33fe4bf7d273228373c88f"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "10k Brain Cells from an E18 Mouse (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "e7b477926c1aa6327bb6da1991b8729d"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "1k Heart Cells from an E18 mouse (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v3/heart_1k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "bad9dd23555d5e228f1875351927e59c"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "10k Heart Cells from an E18 mouse (v3 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "158a95a25ea063b6fa04ecf85fad8ba8"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "10k Mouse E18 Combined Cortex, Hippocampus and Subventricular"
#                         " Zone Cells, Single Indexed"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-single-indexed-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_SI_Neuron_10K/SC3_v3_NextGem_SI_Neuron_10K_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "be561d5e1cede800fb8cbc18df3ef865"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v3"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": (
#                         "10k Mouse E18 Combined Cortex, Hippocampus and Subventricular"
#                         " Zone Cells, Dual Indexed"
#                     )
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/10-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-dual-indexed-3-1-standard-4-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_Neuron_10K/SC3_v3_NextGem_DI_Neuron_10K_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "f9b3896ee55dbcc50e5e1fbc039cdb88"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v2"}},
#             {"key": {"string": "reference"}, "value": {"string": "human2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {"string": "V2 1k PBMCs from a Healthy Donor (v2 chemistry)"},
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-2-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "8dc07a65163a9c1bdfa3c0f1fc7ec653"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v2"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": "V2 1k Brain Cells from an E18 Mouse (v2 chemistry)"
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v2/neuron_1k_v2_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "faae0693e1ccb9806c27f64724be8479"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
#     {
#         "map": [
#             {"key": {"string": "chemistry"}, "value": {"string": "v2"}},
#             {"key": {"string": "reference"}, "value": {"string": "mm10-2020A"}},
#             {
#                 "key": {"string": "dataset_name"},
#                 "value": {
#                     "string": "V2 1k Heart Cells from an E18 mouse (v2 chemistry)"
#                 },
#             },
#             {
#                 "key": {"string": "dataset_url"},
#                 "value": {
#                     "string": "https://www.10xgenomics.com/resources/datasets/1-k-heart-cells-from-an-e-18-mouse-v-2-chemistry-3-standard-3-0-0"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_url"},
#                 "value": {
#                     "string": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v2/heart_1k_v2_fastqs.tar"
#                 },
#             },
#             {
#                 "key": {"string": "fastq_MD5sum"},
#                 "value": {"string": "943b38dbe4c820e5803085ecb7e4752a"},
#             },
#             {"key": {"string": "delete_fastq"}, "value": {"string": "1"}},
#             {"key": {"string": "feature_barcode_csv_url"}, "value": {"string": ""}},
#             {
#                 "key": {"string": "multiplexing_library_csv_url"},
#                 "value": {"string": ""},
#             },
#         ]
#     },
# ]]


try:
    shutil.rmtree(Path(".latch/task-outputs"))
except:
    pass

channel_vals = (
    [
        [{"string": "v3"}],
        [{"string": "human2020A"}],
        [{"string": "500 Human PBMCs, 3' LT v3.1, Chromium Controller"}],
        [{
            "string": "https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-controller-3-1-low-6-1-0"
        }],
        [{
            "string": "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_Controller/500_PBMC_3p_LT_Chromium_Controller_fastqs.tar"
        }],
        [{"string": "5f080c6082f11ea9fc6448482e6fb590"}],
        [{"string": "1"}],
        [{"null": None}],
        [{"null": None}],
        [{"path": "/root/work/31/43ee3c9a970fe1ed8f88ec12b4e32e/3M-february-2018.txt"}],
        [{
            "path": (
                "/root/work/31/43ee3c9a970fe1ed8f88ec12b4e32e/splici_fl85_t2g_3col.tsv"
            )
        }],
        [{
            "path": "/root/work/31/43ee3c9a970fe1ed8f88ec12b4e32e/5f080c6082f11ea9fc6448482e6fb590_alevin_map"
        }],
    ],
)
channel_vals = sum(channel_vals, [])

subprocess.run(
    [
        ".latch/bin/nextflow",
        "run",
        "main.nf",
        "--permitlist",
        "input_files/pl_sheet.tsv",
        "--reference",
        "input_files/ref_sheet.tsv",
        "--sample",
        "input_files/sample_sheet.tsv",
    ],
    env={
        **os.environ,
        "LATCH_INCLUDE_META": (
            '{"path": "./modules/salmon_map.nf", "alias": "salmon_map_rad", "name":'
            ' "salmon_map_rad"}'
        ),
        "LATCH_RETURN": (
            '["{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"chemistry\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"chemistry\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"reference\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"reference\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"dataset_name\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"dataset_name\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"dataset_url\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"dataset_url\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"fastq_url\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"fastq_url\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"fastq_MD5sum\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"fastq_MD5sum\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"delete_fastq\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"delete_fastq\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"feature_barcode_csv_url\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"feature_barcode_csv_url\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"multiplexing_library_csv_url\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"multiplexing_library_csv_url\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"pl_path\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"pl_path\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"t2g_path\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"t2g_path\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"map_dir_path\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"map_dir_path\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"map_rad_path\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"map_rad_path\\"}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"unmapped_file_path\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"salmon_map_rad\\"},\\"property\\":\\"out\\"}},\\"property\\":\\"unmapped_file_path\\"}}}},\\"labels\\":[]}}"]'
        ),
        "LATCH_EXPRESSION": Path("map_expr.json").read_text(),
        "LATCH_PARAM_VALS": json.dumps(channel_vals),
    },
    check=True,
)


out_channels = {}
files = [Path(f) for f in glob.glob(".latch/task-outputs/*.json")]

for file in files:
    out_channels[file.stem] = file.read_text()

print(out_channels)

# paths = []
# for channel in out_channels.values():
#     for param in json.loads(channel):
#         _extract_paths(param, paths)

# for path in paths:
#     print(path)

# print(out_channels)
# upload_files(
#     {k: json.loads(v) for k, v in out_channels.items()},
#     LatchDir("latch://1721.account/Nextflow Outputs/local-testing"),
# )

# return Dataclass_9_post(res=out_channels.get(f"res", ""))
