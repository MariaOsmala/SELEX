
cd /scratch/project_2013895/SELEX/data/


#tar -cjvf Jolma2013.tar.bz Jolma2013/submitted_ftp

#Move to allas

#Creates a new bucket with public access and uploads the data to the bucket. 

#Command a-publish creates the bucket and uploads the selected files into it. 
#Parameter -b is used to define the name for the bucket, in this case TFBS-project-public.

module load allas
allas-conf project_2013895

#Check what there is in allas
a-list -d Jolma2013 | wc -l #2727 DONE!
ls /scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp | wc -l #2726
a-list -d Jolma2015 | wc -l #7049
ls /scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp | wc -l #14847 interactive
a-list -d Morgunova2015 | wc -l #3 DONE!
ls /scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp | wc -l #2
a-list -d Nitta2015 | wc -l #1244 DONE!
ls /scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp | wc -l #1241

#a-delete SELEX-public
#SELEX-public_segments
a-list -d Xie2025 | wc -l #2527 DONE!
ls /scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp | wc -l #2526 
a-list -d Yin2017 | wc -l #5513
ls /scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp | wc -l #7452 #interactive2
a-list -d input_libraries | wc -l #1030 DONE!
ls /scratch/project_2013895/SELEX/data/input_libraries/submitted_ftp | wc -l #1029


#a-publish -b SELEX-public Jolma2013.tar.bz Jolma2013.tar.bz

cd /scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp

for file in *; do
  echo "Processing $file"
  a-publish -b Jolma2013 $file
done

#Processing ERF_E2F4_1_AAC_TTACAG40NGTG.fastq.gz
cd /scratch/project_2013895/SELEX/data/Jolma2015/submitted_ftp

files=(*)
echo ${files[7049]} #HOXC10_CREB3L1_3_AZ_TGCTTG40NTGC.fastq.gz

mapfile -t included_files < /projappl/project_2013895/SELEX/data-to-allas/Jolma2015.txt

# Turn included_files into an associative array for fast lookup
declare -A incmap
for f in "${included_files[@]}"; do
    incmap["$f"]=1
done

# Filter files, keep original order
filtered=()
for f in "${files[@]}"; do
    if [[ ${incmap["$f"]+yes} ]]; then
        filtered+=("$f")
    fi
done

# Result
printf '%s\n' "${filtered[@]}"

#What there already is in Allas
mapfile -t added < <(a-list -d Jolma2015)

# build lookup table from added
declare -A in_added
for f in "${added[@]}"; do
    in_added["$f"]=1
done

# collect those in filtered but not in added
not_in_added=()
for f in "${filtered[@]}"; do
    if [[ -z ${in_added["$f"]+yes} ]]; then
        not_in_added+=("$f")
    fi
done

# print result
printf '%s\n' "${not_in_added[@]}"




#for file in *; do

#ls | awk '/ERF_E2F4_1_AAC_TTACAG40NGTG.fastq.gz/ {found=1} found'
files=$(ls | awk '/ERF_E2F4_1_AAC_TTACAG40NGTG.fastq.gz/ {found=1} found')

#ls | awk '/FOXO1_SOX17_1_AS_TGACTA40NGGC.fastq.gz/ {found=1} found'
files=$(ls | awk '/FOXO1_SOX17_1_AS_TGACTA40NGGC.fastq.gz/ {found=1} found')

files=$(ls | awk '/HOXA3_PITX1_3_AY_TACTAG40NGGA.fastq.gz/ {found=1} found')


for file in $files; do
  echo "Processing $file"
  a-publish -b Jolma2015 $file
done

for file in ${not_in_added[@]}; do
  echo "Processing $file"
  a-publish -b Jolma2015 $file
done


#https://a3s.fi/Jolma2015/293FT_MEIS1_sorted.bam.fastq.gz

cd /scratch/project_2013895/SELEX/data/input_libraries/submitted_ftp

#ls | awk '/ZeroCycle_TTTGTT40NTTAG_0_0.fastq.gz/ {found=1} found'
#This is finished
for file in *; do
  echo "Processing $file"
  a-publish -b input_libraries $file
done

#This is finished
cd /scratch/project_2013895/SELEX/data/Morgunova2015/submitted_ftp
for file in *; do
  echo "Processing $file"
  a-publish -b Morgunova2015 $file
done

#ls | awk '/KY_TTTGTT40NTAT_4.fastq.gz/ {found=1} found'
#This is finished
cd /scratch/project_2013895/SELEX/data/Nitta2015/submitted_ftp
for file in *; do
  echo "Processing $file"
  a-publish -b Nitta2015 $file
done


cd /scratch/project_2013895/SELEX/data/Xie2025/submitted_ftp
for file in *; do
  echo "Processing $file"
  a-publish -b Xie2025 $file
done



cd /scratch/project_2013895/SELEX/data/Yin2017/submitted_ftp

files=(*)


mapfile -t included_files < /projappl/project_2013895/SELEX/data-to-allas/Yin2017.txt

# Turn included_files into an associative array for fast lookup
declare -A incmap
for f in "${included_files[@]}"; do
    incmap["$f"]=1
done

# Filter files, keep original order
filtered=()
for f in "${files[@]}"; do
    if [[ ${incmap["$f"]+yes} ]]; then
        filtered+=("$f")
    fi
done

# Result
printf '%s\n' "${filtered[@]}"

#What there already is in Allas
mapfile -t added < <(a-list -d Yin2017)

# build lookup table from added
declare -A in_added
for f in "${added[@]}"; do
    in_added["$f"]=1
done

# collect those in filtered but not in added
not_in_added=()
for f in "${filtered[@]}"; do
    if [[ -z ${in_added["$f"]+yes} ]]; then
        not_in_added+=("$f")
    fi
done

# print result
printf '%s\n' "${not_in_added[@]}"






#ls | awk '/LHX8_FL_3_KW_TTGCGA40NGTA.fastq.gz/ {found=1} found'
files=$(ls | awk '/LHX8_FL_4_KW_TAGCTG40NCTA.fastq.gz/ {found=1} found')
files=$(ls | awk '/PRRX2_FL_3_KV_TGCCAA40NTAG.fastq.gz/ {found=1} found')

for file in $files; do
#for file in *; do
  echo "Processing $file"
  a-publish -b Yin2017 $file
done


for file in ${not_in_added[@]}; do
  echo "Processing $file"
  a-publish -b Yin2017 $file
done



https://a3s.fi/Jolma2013/ARNTL_TCAAAA20NCG_W_1.fastq.gz




tmp=$(ls | head)

a-publish -b Jolma2013 NOTO_TGCGTT30NTGC_AI_1.fastq.gz

a-publish -b Jolma2013 --input-list Alx1_TAAAGC20NCG_Z_1.fastq.gz Alx1_TAAAGC20NCG_Z_2.fastq.gz Alx1_TAAAGC20NCG_Z_3.fastq.gz Alx1_TAAAGC20NCG_Z_4.fastq.gz 

a-publish -b SELEX-public --input-list $tmp


#a-publish copies a file to Allas into a bucket that can be publicly accessed. 
#Thus, anyone with the address (URL) of the  uploaded data object can read and download the data with a web browser or tools like wget and curl. 
#a-publish works mostly like a-put but there are some differences: 

#1) a-publish can upload only files, not directories. 
#2) files are not compressed but they uploaded as they are. 
#3) the access control of the target bucket is set so that it is available in read-only mode to the internet.

#The basic syntax of the command is:

#a-publish file_name

#By default, the file is uploaded to a bucket username-projectNumber-pub. You can define other bucket names too using option -b
#but you should note that this command will make all data in the bucket publicly accessible, 
#including data that has been previously uploaded to the bucket.

#The public URL to a data object is:

#https://a3s.fi/username-projectNumber-pub/object_name

https://a3s.fi/SELEX-public/Jolma2013/NOTO_TGCGTT30NTGC_AI_1.fastq.gz

#An object uploaded with a-publish can be removed from Allas with command a-delete.

#a-publish options:

# -b, --bucket       Use the defined bucket instead of the default bucket name
# -o, --os_file      Define alternative name for the object that will be created  

# -i, --index        (static/dynamic).  By default a-publish creates a static index 
#                   file that includes the objects that are in the target bucket when 
#                    the command is executed.
#                     With setting --index dynamic the command adds a JavaScript based 
#                    index file to the bucket. With this option the index.html page 
#                    lists the objects that are available in the bucket in the time when 
#                    this page is accessed. This dynamic indexing tool can list
#                    only up to 1000 files.
                    
                    
 #--input-list       List of files to be uploaded.    
 


# ps
# Jos teet noista tuloksista bed muotoisen,  
# pakkaat bgzip:illä, indeksoit tabix:illa, lataat paketin ja indeksin altaalle,  
# jaat public:ille, ja ajat alla olevan komennon (sopivasti komentoriviä muokaten,
# s3cors_policy.xml tiedosto liittenä),  niin tiedostoja voi käyttää kätevästi igv:llä (Load from URL..) 
# indeksoituna ilman että tarvii ladata kaikkea.

# s3cmd setcors /home/kpalin/src/myscripts/s3cors_policy.xml  s3://240704_ont_pipe_chm13v2_Gp5d_ELF2_Hia5_1