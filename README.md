# scmultiome_preprocessing
The automated single cell RNA and single cell ATAC preprocessing pipeline using Seurat and Signac.
This will also run ChromVar for Transcription factor enrichment

## To Run the preprocessing of scmultiome
      
      source("scmultiome_Seurat_processing.R")
      savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/final_seurat_obj/naive_resting/"
      
      rest_act <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/final_seurat_obj/final_obj_res_act_old_young.RDS") # object should have Age
      rest_act@meta.data$Age <- gsub("#.*|#.*","",rownames(rest_act@meta.data))
      
      cellnames <- read.delim("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scmultiome/Huimin_and_DrGoronzy/resting_cellnames.txt")
      
      scmultiome_pre(obj = rest_act,savedir = savedir, objname = "naive_resting", subset_cellnames = cellnames[,1])

![Screenshot 2024-04-01 at 6 19 45 PM](https://github.com/Ajaingithub/scmultiome_preprocessing/assets/37553954/fa15b205-4e02-4c0f-88c7-a2637071deef)

![Screenshot 2024-04-01 at 6 21 24 PM](https://github.com/Ajaingithub/scmultiome_preprocessing/assets/37553954/176034a6-01d2-488e-9347-1fc8ef92829c)
