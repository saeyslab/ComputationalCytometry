install.packages("devtools")
install.packages("Rtsne")
install.packages("tidyverse")
install.packages("openxlsx")
install.packages("mvtnorm")
install.packages("reshape2")
install.packages("gridExtra")
install.packages("R.utils")
install.packages("IsolationForest", repos="http://R-Forge.R-project.org")

# Bioconductor packages , to be installed through BiocManager:

install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("flowAI")
BiocManager::install("flowDensity")
BiocManager::install("CytoML")
BiocManager::install("flowStats")

# FlowSOM and flowSOM_workshop packages, to be installed directly from the
# github repository:

devtools::install_github("SofieVG/FlowSOM", ref = "devel")
# If you're on an older bioconductor version, try 
#devtools::install_github("saeyslab/FlowSOM") 

devtools::install_github("saeyslab/FlowSOM_workshop")
devtools::install_github("saeyslab/CytoNorm")
