cd /projappl/project_2013895/softwares
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz

tar -xvzf jellyfish-2.3.1.tar.gz

module load gcc
cd jellyfish-2.3.1
./configure 
make -j 4
make install

PATH=/projappl/project_2013895/softwares/jellyfish-2.3.1/bin:$PATH