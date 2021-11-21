pip install gdown
cd
#gdown https://drive.google.com/uc?id=1nIy0DZl1p_nR-bAvYtfkiI_-B68Q7Xzy
wget "https://myfiles.lshtm.ac.uk/rest/files/public/links/8a8c80b67d24d440017d3dea13b551c7?passKey=-3471261619226576060&shareId=9049" -O course_data.tar.gz
tar -xzvf course_data.tar.gz
rm course_data.tar.gz
mv course_data data
