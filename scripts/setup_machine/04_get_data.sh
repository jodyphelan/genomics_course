pip install gdown
cd
#gdown https://drive.google.com/uc?id=1nIy0DZl1p_nR-bAvYtfkiI_-B68Q7Xzy
# wget "https://myfiles.lshtm.ac.uk/filr/public-link/file-download/8a8c80b78443190b0184c3dc68e77a06/11340/6751946147449931689/course_data.tar.gz" -O course_data.tar.gz
wget https://genomics.lshtm.ac.uk/course_data_v1.2.tgz
tar -xzvf course_data_v1.2.tgz
rm course_data.tar.gz
cd ~/data/nanopore_activity/kraken
# wget https://myfiles.lshtm.ac.uk/filr/public-link/file-download/8a8c80b78443190b0184cf6e37d4346a/11347/-4475221956266795664/taxdmp.zip
# unzip taxdmp.zip
# rm taxdmp.zip
