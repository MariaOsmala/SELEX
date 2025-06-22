
cd /scratch/project_2013895/SELEX/data/


tar -cjvf Jolma2013.tar.bz Jolma2013/submitted_ftp

#Move to allas

#Creates a new bucket with public access and uploads the data to the bucket. 

#Command a-publish creates the bucket and uploads the selected files into it. 
#Parameter -b is used to define the name for the bucket, in this case TFBS-project-public.

module load allas
allas-conf project_2013895


a-publish -b SELEX-public Jolma2013.tar.bz Jolma2013.tar.bz

cd /scratch/project_2013895/SELEX/data/Jolma2013/submitted_ftp

for file in *; do
  echo "Processing $file"
  a-publish -b Jolma2013 $file
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