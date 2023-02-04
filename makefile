.PHONY: install

install:
	
	@ # https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox#4-explore-virtual-box
	# sudo apt update && sudo apt upgrade -y
	# sudo snap refresh
	
	@ # install R
	@ # https://stackoverflow.com/a/38698271
	@ # https://stackoverflow.com/questions/10608004/auto-install-packages-from-inside-makefile
	[ -z "`dpkg -l | grep libcurl4-openssl-dev`" ] && sudo apt-get install libcurl4-openssl-dev || echo
	[ -z "`dpkg -l | grep libssl-dev`" ] && sudo apt-get install libssl-dev || echo
	[ -z "`dpkg -l | grep libfontconfig1-dev`" ] && sudo apt-get install libfontconfig1-dev || echo
	[ -z "`dpkg -l | grep libxml2-dev`" ] && sudo apt-get install libxml2-dev || echo
	[ -z "`dpkg -l | grep libharfbuzz-dev`" ] && sudo apt-get install libharfbuzz-dev || echo
	[ -z "`dpkg -l | grep libfribidi-dev`" ] && sudo apt-get install libfribidi-dev || echo
	[ -z "`dpkg -l | grep libfreetype6-dev`" ] && sudo apt-get install libfreetype6-dev || echo
	[ -z "`dpkg -l | grep libpng-dev`" ] && sudo apt-get install libpng-dev || echo
	[ -z "`dpkg -l | grep libtiff5-dev`" ] && sudo apt-get install libtiff5-dev || echo 
	[ -z "`dpkg -l | grep libjpeg-dev`" ] && sudo apt-get install libjpeg-dev || echo
	[ -z "`dpkg -l | grep r-base-dev`" ] && sudo apt-get install r-base-dev || echo
	
	@ # install dependencies
	@ # https://stackoverflow.com/questions/6907937/how-to-install-dependencies-when-using-r-cmd-install-to-install-r-packages
	sudo Rscript -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
	sudo Rscript -e 'devtools::install_deps(".")'

	@ # build package
	sudo R CMD INSTALL .

	
	

