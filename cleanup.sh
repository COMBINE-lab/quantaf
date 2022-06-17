#!/bin/bash

./nextflow main.nf -profile test
rm -rf work
rm -rf nf_pipeline
rm .nextflow.log*
rm pipeline_trace.txt*