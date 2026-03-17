cd /scratch/project_2013895/SELEX 

tar -zcvf sequence-to-affinity-data.tar.gz sequence-to-affinity-data/ #16G -> 982M

module load allas
allas-conf --mode S3

s3cmd mb s3://SELEX/

#List all buckets in a project:

s3cmd ls

#Upload a file to a bucket:

s3cmd put sequence-to-affinity-data.tar.gz s3://SELEX/sequence-to-affinity-data.tar.gz

#List all objects in a bucket
s3cmd ls s3://SELEX/

#Make the bucket public
s3cmd setacl --acl-public s3://SELEX

#Make file publice
s3cmd setacl --acl-public s3://SELEX/sequence-to-affinity-data.tar.gz

#The syntax of the URL of the file:

https://a3s.fi/SELEX/sequence-to-affinity-data.tar.gz

#Display information about a bucket:
s3cmd info s3://SELEX

#Display information about an object:
s3cmd info s3://SELEX/sequence-to-affinity-data.tar.gz

#Check your project's object storage usage:
s3cmd du -H

#Download an object:
s3cmd get s3://my_bucket/my_file new_file_name

#Download an entire bucket:
s3cmd get -r s3://my_bucket/

#Copy an object to another bucket. Note that should use these commands only 
# for objects that were uploaded to Allas with S3 protocol:
s3cmd cp s3://sourcebucket/objectname s3://destinationbucket

s3cmd cp s3://bigbucket/bigfish s3://my-new-bucket
#remote copy: 's3://bigbucket/bigfish' -> 's3://my-new-bucket/bigfish'

#Rename the file while copying it:
s3cmd cp s3://bigbucket/bigfish s3://my-new-bucket/newname
#remote copy: 's3://bigbucket/bigfish' -> 's3://my-new-bucket/newname'

#Delete an object:
s3cmd del s3://my_bucket/my_file
#Delete a bucket:
s3cmd rb s3://my_bucket
#Note: You can only delete empty buckets.

#s3cmd and public objects
#In this example, the object salmon.jpg in the pseudo folder fishes is made public:

#$ s3cmd put fishes/salmon.jpg s3://my_fishbucket/fishes/salmon.jpg -P
#Public URL of the object is: https://a3s.fi/my_fishbucket/fishes/salmon.jpg


#using the command s3cmd setacl, you can make the file publicly available.

#First make the fish bucket public:

s3cmd setacl --acl-public s3://fish-bucket
Then make the zebrafish genome file public:


s3cmd setacl --acl-public s3://fish-bucket/zebrafish.tgz
The syntax of the URL of the file:

s3cmd setacl --acl-public s3://random_pe_manu_preprocess_1911/
s3cmd setacl --acl-public s3://random_pe_manu_preprocess_1911/fastq/InputRandomPE/InputRandomPE_R1.fastq.gz
s3cmd setacl --acl-public s3://random_pe_manu_preprocess_1911/fastq/InputRandomPE/InputRandomPE_R2.fastq.gz
