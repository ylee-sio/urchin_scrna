setwd("~/projects")
library(tidyverse)
library(gtools)

list.files()
selection = as.numeric(ask("Pick project: "))
project_selection = list.files()[selection]
setwd(project_selection)

print(paste0("Set current working directory to ", project_selection))
list.files()
selection = as.numeric(ask("Source script: "))
source(list.files()[selection])

t.out = "temp_outputs"

