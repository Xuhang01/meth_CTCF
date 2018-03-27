### run in root folder:  bash ./utils/preprocessing/splitting_meth_file.sh  OV-US
proj_code=$1
input_file=/f/hang/PhD_Project/Data/Database/ICGC/ICGC_Project/$proj_code/meth_array.$proj_code.tsv.gz
# add file name column (int(cg********[2:])/100000 one file for each 100000 interval)
output_file=/f/hang/PhD_Project/Project/meth_anchor/data/split_raw/raw_add_columns/meth_array.$proj_code.tsv
time gunzip -c $input_file | parallel --pipe -k -q awk -F $'\t' -v proj_code="$proj_code" \
'{if(($7=="HumanMethylation450_after_2011_08_02") && match($8, /^cg/)){v=sprintf("%d",substr($8,3)/100000);print proj_code"_"v"\t"$0}}' \
 > $output_file
# split file to different files
cd /f/hang/PhD_Project/Project/meth_anchor/data/split_raw/
mkdir $proj_code
cd $proj_code
time cat $output_file  | awk  -F $'\t' '{print >$1".tsv"}'
ls $proj_code\_*.tsv | parallel gzip
# gzip the add_column file
gzip $output_file


### run in root folder:  /f/hang/PhD_Project/Project/meth_anchor/
# nohup bash ./utils/preprocessing/splitting_meth_file.sh BLCA-US > ./data/split_raw/BLCA-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh BRCA-US > ./data/split_raw/BRCA-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh CESC-US > ./data/split_raw/CESC-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh COAD-US > ./data/split_raw/COAD-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh GBM-US  > ./data/split_raw/GBM-US.nohup.txt   &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh KIRC-US > ./data/split_raw/KIRC-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh KIRP-US > ./data/split_raw/KIRP-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh LAML-US > ./data/split_raw/LAML-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh LIHC-US > ./data/split_raw/LIHC-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh LUAD-US > ./data/split_raw/LUAD-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh LUSC-US > ./data/split_raw/LUSC-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh OV-US   > ./data/split_raw/OV-US.nohup.txt    &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh PAAD-US > ./data/split_raw/PAAD-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh PRAD-US > ./data/split_raw/PRAD-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh READ-US > ./data/split_raw/READ-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh SKCM-US > ./data/split_raw/SKCM-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh STAD-US > ./data/split_raw/STAD-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh THCA-US > ./data/split_raw/THCA-US.nohup.txt  &
# nohup bash ./utils/preprocessing/splitting_meth_file.sh UCEC-US > ./data/split_raw/UCEC-US.nohup.txt  &