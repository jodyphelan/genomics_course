mkdir ~/git
cd ~/git
git clone https://github.com/julibeg/dotfiles.git
find ~/git/dotfiles/ -maxdepth 1 -type f -exec ln -s {} ~/ \;
git clone https://github.com/jodyphelan/genomics_course.git
#source ~/git/genomics_course/scripts/setup_machine/*sh
