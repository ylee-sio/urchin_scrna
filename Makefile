miniconda-env:
	mkdir ~/projects
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh -P ~/projects
	bash  ~/projects/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh
	source ~/.bashrc
	source ~/.bash_profile
	conda create --name rconda-seurat
	conda activate rconda-seurat
	conda install r-base r-tidyverse r-seurat
	
datasets:
	mkdir ~/projects/urchin_scrna/data_sources
	wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149221/suppl/GSE149221_SpInteg.rds.gz -P  ~/projects/urchin_scrna/data_sources
	wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155427/suppl/GSE155427_Sp48MO.rds.gz -P  ~/projects/urchin_scrna/data_sources
	wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155427/suppl/GSE155427_Sp48and72.rds.gz -P  ~/projects/urchin_scrna/data_sources
	gunzip  ~/projects/urchin_scrna/data_sources/*gz

 
