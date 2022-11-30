mkdir -p ~/.vim/{backup,swap,undo}
mkdir ~/git
cd ~/git
#git clone https://github.com/julibeg/dotfiles.git
#find ~/git/dotfiles/ -maxdepth 1 -type f -exec ln -s {} ~/ \;
git clone https://github.com/jodyphelan/genomics_course.git
source <(cat ~/git/genomics_course/scripts/setup_machine/0*sh)
ln -s ~/git/genomics_course/scripts/update_course.sh ~/bin
chmod +x ~/bin/update_course.sh
