#!/bin/bash

# docker run -itv .:/workspace -v /data/GL/database:/database  namxle/scip-backend:2.0.0 bash
# cd /workspace/scip && rm -rf * && /bin/bash /workspace/run.sh 3959

scip_workspace_dir="/scip-workspace"

sample_id=${1}
dir="/workspace/${sample_id}"

cnv_vcf="${dir}/${sample_id}.cnv.vcf"
deduped_bam="${dir}/${sample_id}.deduped.bam"
deduped_bam_bai="${dir}/${sample_id}.deduped.bam.bai"

genome_name=hg19

scip_sample_file="sample_id.txt"
scip_sample_id="${sample_id}-001-001"

scip_unfiltered_cnv="${scip_sample_id}.unfiltered_CNV.txt"

scip_alignment_dir="${scip_workspace_dir}/alignment/${sample_id}"
scip_input_bam="${sample_id}.bam"
scip_input_bam_bai="${sample_id}.bam.bai"

scip_user_data_dir="user_data"
scip_app_tmp_f_dir="app_temp_file"

scip_summary="${scip_sample_id}.${genome_name}.pipeline_summary.txt"

scip_pl_script="SCIP_backend_${genome_name}.pl"

scip_user_data_output="${scip_user_data_dir}/"
scip_app_temp_file_output="${scip_app_tmp_f_dir}/"
scip_interface_config="interface_config.txt"
scip_log="scip.log"

# Original working directory
OWD=$PWD

# Generate some directories
mkdir -p ${scip_user_data_dir} ${scip_app_tmp_f_dir}

mkdir -p ${scip_workspace_dir}/${scip_user_data_dir}
mkdir -p ${scip_workspace_dir}/${scip_app_tmp_f_dir}

# Prepare CNV file
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ALT\n' ${cnv_vcf} |
    awk -F"\t" 'BEGIN{OFS="\t"}{sub(/^chr/,"",$1)}1' |
    awk -F"\t" -v scip_sample_id=${scip_sample_id} '{ \
        split($4, a, ">"); \
        split(a[1], b, "<"); \
        print $1"\t"$2"\t"$3"\t"b[2]"\t.\t.\t"scip_sample_id; \
    }' \
    > ${scip_unfiltered_cnv}

# ls -la ${scip_workspace_dir} >> test.log

cp $OWD/${scip_unfiltered_cnv} ${scip_workspace_dir}/${scip_user_data_dir}/${scip_unfiltered_cnv}

# Prepare alignment data
mkdir -p ${scip_alignment_dir}
ln -s ${deduped_bam} ${scip_alignment_dir}/${scip_input_bam}
ln -s ${deduped_bam_bai} ${scip_alignment_dir}/${scip_input_bam_bai}

# Prepare sample ID data
echo "${sample_id}	${scip_sample_id}" > ${scip_sample_file}
cp ${scip_sample_file} ${scip_workspace_dir}/${scip_sample_file}

# Run SCIP backend
cd ${scip_workspace_dir}
perl ${scip_pl_script} -n ${scip_sample_id} -@ 40 > ${scip_log}

cp -r ${scip_workspace_dir}/${scip_user_data_dir} $OWD
cp -r ${scip_workspace_dir}/${scip_app_tmp_f_dir} $OWD
cp ${scip_workspace_dir}/${scip_log} $OWD

# Check the result of SCIP
#if [[ $(wc -l < "${scip_workspace_dir}/${scip_user_data_dir}/${scip_summary}") -eq 0 ]]; then
#    echo "File ${scip_summary} has no lines. Exiting with status 1."
#    exit 1
#fi

# Create interface config file
cd $OWD
echo "LIST_NAME	${scip_sample_id}.${genome_name}" > ${scip_interface_config}
echo "TEMP_FILE_DIR	/highspeed-data/samples/${sample_id}/SCIP/app_temp_file/" >> ${scip_interface_config}
echo "ROOT_DIR	/highspeed-data/samples/${sample_id}/SCIP/" >> ${scip_interface_config}
echo "USER	test" >> ${scip_interface_config}