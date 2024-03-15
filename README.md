# scmultiome_preprocessing
The automated single cell RNA and single cell ATAC preprocessing pipeline using Seurat and Signac.

## To Run the preprocessing of scmultiome
      
      source("scmultiome_Seurat_processing.R")
      savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/final_seurat_obj/naive_resting/"
      
      rest_act <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/final_seurat_obj/final_obj_res_act_old_young.RDS") # object should have Age
      rest_act@meta.data$Age <- gsub("#.*|#.*","",rownames(rest_act@meta.data))
      
      cellnames <- read.delim("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/resting_cellnames.txt")
      
      scmultiome_pre(obj = rest_act,savedir = savedir, objname = "naive_resting", subset_cellnames = cellnames[,1])
