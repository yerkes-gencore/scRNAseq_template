##############################################
# Starting from freshly forked template repo #
##############################################

## Run this script when first initializing a new project from the template
## Or if scrapping a corrupt renv setup, after running `renv::deactivate(clean=TRUE)`

# Make sure R package manager repo is up-to-date
# options(repos = "https://packagemanager.rstudio.com/all/__linux__/focal/latest")
renv::install("remotes")
## If using R 4.3.1 [[ NOTE: You can set rebuild = FALSE (Default) since I've done this and it's in the global cache ]]
remotes::install_version("Matrix", version = "1.6-4") # see https://github.com/plger/scDblFinder/issues/91#issuecomment-1825482069
renv::install("irlba", rebuild = TRUE, prompt = FALSE)
renv::install("Seurat@5.0.3", rebuild = TRUE, prompt = FALSE) # see https://github.com/satijalab/seurat/issues/8916#issuecomment-2164844005
renv::install("SeuratObject@5.0.1", rebuild = TRUE, prompt = FALSE) # also see https://github.com/satijalab/seurat/issues/8916#issuecomment-2164844005
## Error: package or namespace load failed for ‘celldex’ in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]): namespace ‘DBI’ 1.1.3 is already loaded, but >= 1.2.0 is required
renv::install("DBI@1.2.0")

## Install the gencore packages.
renv::install("yerkes-gencore/gencoreSC")
renv::install("yerkes-gencore/gencoreBulk")

## Install loupeR from github repo
renv::install("10XGenomics/loupeR")

## 
renv::install(c("boot", "nnet", "writexl"))

## Initialize, and setup renv to work with bioconductor
renv::init(bioconductor = TRUE)

## Seurat v5: If using Seurat v5 you should install the following suggested packages as well
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
renv::install(c("BPCells", "presto", "glmGamPoi", "hdf5r"))

## Run renv::status() and resolve any issues before snapshotting
renv::status()

## Snapshot the package versions so that they can be reproduced
renv::snapshot()

# #############################################################
# # After cloning an existing project into a new git "remote" #
# #############################################################
#
# ## Initialize a project with the packages from the other remote (if the lockfile was up-to-date prior to cloning the repo)
# renv::restore()
#
# ## Run renv::status() and resolve any issues before snapshotting
# renv::status()
#
# ## Snapshot the package versions so that they can be reproduced
# renv::snapshot()
