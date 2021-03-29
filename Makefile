miniconda-env:
	mkdir ~/projects
	mkdir ~/plots
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh -P ~/projects
	bash  ~/projects/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh
	source ~/.bashrc
	source ~/.bash_profile
	conda create --name rconda-seurat
	conda activate rconda-seurat
	conda install r-base r-tidyverse r-seurat parallel
	
datasets:
	mkdir ~/projects/urchin_scrna/data_sources
	cat installation/data_download_links.txt | parallel -j 8 wget -P ~/projects/urchin_scrna/data_sources/
	gunzip  ~/projects/urchin_scrna/data_sources/*

plotset:
	mkdir ~projects/urchin_scrna/plots/plotset{#} {}
 
