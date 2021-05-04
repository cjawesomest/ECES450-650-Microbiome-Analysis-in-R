## Installing packages
.cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages, version = "3.12") #Different from the website's "3.9" since I'm using R 4.0.5
devtools::install_github("adw96/breakaway")
devtools::install_github(repo = "UVic-omics/selbal")

## Loading Libraries
library(tidyverse); packageVersion("tidyverse")    
library(phyloseq); packageVersion("phyloseq")   
library(DESeq2); packageVersion("DESeq2")
library(microbiome); packageVersion("microbiome") 
library(vegan); packageVersion("vegan")  
library(picante); packageVersion("picante") 
library(ALDEx2); packageVersion("ALDEx2") 
library(metagenomeSeq); packageVersion("metagenomeSeq") 
library(HMP); packageVersion("HMP")  
library(dendextend); packageVersion("dendextend")  
library(selbal); packageVersion("selbal")  
library(rms); packageVersion("rms")
library(breakaway); packageVersion("breakaway") 