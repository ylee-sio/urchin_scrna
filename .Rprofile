setwd("~/projects")
library(parallel)
library(tidyverse)
library(gtools)
library(Seurat)

# displays options
list.files()

# sets project working directory
selection = as.numeric(ask("Pick project: "))
project_selection = list.files()[selection]
setwd(project_selection)
print(paste0("Set current working directory to ", project_selection))
print("************************")

source("r/source.R")
system.time(source(paste0("r/",source_dir[selection])))

