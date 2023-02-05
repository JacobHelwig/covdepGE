# to install covdepGE from CRAN, run install.packages("covdepGE", repos="https://cloud.r-project.org")

.PHONY: install

install:

	@ # https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#4-explore-virtual-box
	sudo apt update && sudo apt upgrade -y
	sudo snap refresh

	@ # install R
	@ # https://stackoverflow.com/a/38698271
	sudo apt-get install libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev r-base-dev

	@ # install dependencies
	@ # https://stackoverflow.com/questions/6907937/how-to-install-dependencies-when-using-r-cmd-install-to-install-r-packages
	sudo R -vanilla -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
	sudo R -vanilla -e 'devtools::install_deps(".")'

	@ # build package
	sudo R CMD INSTALL .




