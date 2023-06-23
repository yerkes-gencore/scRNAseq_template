install.packages(c('remotes'))
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remotes::install_github('yerkes-gencore/gencoreSC')

install.packages('BiocManager')
BiocManager::install("BiocParallel")
BiocManager::install("MAST")
BiocManager::install("SingleR")

remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github('satijalab/seurat-data')
## need seurat disk and seurat data prior to azimuth
remotes::install_github('satijalab/azimuth')