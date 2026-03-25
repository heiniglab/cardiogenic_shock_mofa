sessionInfo()
.libPaths("/ictstr01/home/icb/bhavishya.nelakuditi/miniconda3/envs/my_jupyter_env/lib/R/library/")
path_to_conda_env = '/ictstr01/home/icb/bhavishya.nelakuditi/miniconda3/envs/my_jupyter_env/lib/R/library/'
print(path_to_conda_env)
library(tidyr, lib.loc = path_to_conda_env)
library(ggplot2, lib.loc = path_to_conda_env)
library(tidyverse, lib.loc = path_to_conda_env)

library(tzdb, lib.loc = path_to_conda_env)
library(readr, lib.loc = path_to_conda_env)
library(backports, lib.loc = path_to_conda_env)
library(ggplot2, lib.loc = path_to_conda_env)
library(ggraph, lib.loc = path_to_conda_env)
library(ggpubr, lib.loc = path_to_conda_env)
library(corrplot, lib.loc = path_to_conda_env)

library(data.table, lib.loc = path_to_conda_env)
library(stringr, quietly = TRUE, verbose = FALSE, lib.loc = path_to_conda_env)
library(dplyr, lib.loc = path_to_conda_env)
library(tidyverse, lib.loc = path_to_conda_env)
library('reshape2', lib.loc = path_to_conda_env)
library(gridExtra, lib.loc = path_to_conda_env)
library(scales, lib.loc = path_to_conda_env)


Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/my_jupyter_env/bin/python")
library(reticulate, quietly = TRUE, verbose = FALSE, lib.loc = path_to_conda_env)
reticulate::use_python("~/miniconda3/envs/my_jupyter_env/bin/python", required=TRUE) 
reticulate::use_condaenv("~/miniconda3/envs/my_jupyter_env", required=TRUE)
reticulate::py_config()

library(ggokabeito, lib.loc = path_to_conda_env)

library(MOFA2, lib.loc = path_to_conda_env)
library(MOFAdata, lib.loc = path_to_conda_env)



