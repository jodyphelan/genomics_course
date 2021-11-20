apt update
apt install -y default-jre xclip
useradd -m -s /bin/bash user
echo user:8eyMByVWFjgy6Na | chpasswd
sudo -i -u user bash << EOF
wget https://raw.githubusercontent.com/jodyphelan/genomics_course/master/scripts/setup_machine/as_user.sh
bash as_user.sh
EOF
