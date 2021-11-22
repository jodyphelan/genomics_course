#!/bin/bash

pw=password

# get packages, desktop, VNC server etc.
apt update && apt upgrade && apt autoremove
apt install -y tigervnc-standalone-server tigervnc-xorg-extension
#apt install -y tigervnc-viewer
DEBIAN_FRONTEND=noninteractive apt install -yq lubuntu-desktop thunar
apt remove -y xscreensaver xfwm4 mutter

# add user
useradd -m user -s /bin/bash
echo user:$pw | chpasswd
echo "shopt -s autocd" >> /home/user/.bashrc

# set up VNC
mkdir /home/user/.vnc
echo $pw  |  vncpasswd -f > /home/user/.vnc/passwd
chown -R user:user /home/user/.vnc
chmod 0600 /home/user/.vnc/passwd

echo "#!/bin/sh" > /home/user/.vnc/xstartup
echo "xrdb $HOME/.Xresources" >> /home/user/.vnc/xstartup
echo "startlxqt &" >> /home/user/.vnc/xstartup
chmod +x /home/user/.vnc/xstartup
sudo -i -u user bash -c 'vncserver -localhost no'

# configure firewall
ufw allow OpenSSH
ufw allow from any to any port 5901
ufw --force enable

# finally prepare machine for the course
apt install -y default-jre xclip
sudo -i -u user bash << EOF
wget https://raw.githubusercontent.com/jodyphelan/genomics_course/master/scripts/setup_machine/as_user.sh
bash as_user.sh
EOF

echo "all done" > /home/user/done
