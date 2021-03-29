basic-structure:
        mkdir ~/projects
        mkdir ~/projects/plots

miniconda-env:
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh -P
	bash  ~/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh
	source ~/.bashrc
	source ~/.bash_profile
	conda create --name rconda-seurat
	conda activate rconda-seurat
	conda install r-base r-tidyverse r-seurat parallel
	
datasets:
	echo "This is not done until you see, 'DONE'."
	mkdir ~/projects/urchin_scrna/data_sources
	cat installation/data_download_links.txt | parallel -j 3 wget -P ~/projects/urchin_scrna/data_sources/
	gunzip  ~/projects/urchin_scrna/data_sources/*
	echo DONE
 
