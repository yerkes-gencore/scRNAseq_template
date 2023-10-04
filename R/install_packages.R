## WIP, idea is to possibly have a script that checks all the packages
## needed for the 'default' analysis, and installs them if needed. This could
## facilitate automation later on

cran_packages <- c(
  'renv',
  'BiocManager',
  'here',
  'tidyverse',
  'Seurat',
  'ggh4x' ## for threshold filtering plots
  #'ggpubr' ## needed?
)

bioc_packages <- c(
  'SingleR',
  'BiocParallel',             ## allow parallelization
  'MAST'                      ## DGE algorithm
)

github_packages <- c(
  'yerkes-gencore/gencoreSC'
  # "mojaveazure/seurat-disk", ## Required for azimuth
  # 'satijalab/seurat-data',   ## Required for azimuth
  # 'satijalab/azimuth'
)

renv::init(restart = FALSE)

for (package in cran_packages) {
  tryCatch(
    {
      if (!require(package, quietly = TRUE, character.only = TRUE)) {
        install.packages(package, character.only = TRUE, clean = TRUE)
        message(paste0(package, ' successfully installed'))
      } else {
        message(paste0(package, ' already installed'))
      }
    },
    error = function(e) {
      message(paste0('Error installing package ', package, '\nOriginal message:'))
      message(e)
    }#,
    # warning=function(w){
    #   message(paste0('Warning installing package ', lib, '\nOriginal message:'))
    #   message(w)
    # }
  )
}

for (package in bioc_packages) {
  tryCatch(
    {
      if (!require(package, quietly = TRUE, character.only = TRUE)) {
        BiocManager::install(package)
        message(paste0(package, ' successfully installed'))
      } else {
        message(paste0(package, ' already installed'))
      }
    },
    error = function(e){
      message(paste0('Error installing package ', package, '\nOriginal message:'))
      message(e)
    }#,
    # warning=function(w){
    #   message(paste0('Warning installing package ', lib, '\nOriginal message:'))
    #   message(w)
    # }
  )
}

for (package in github_packages) {
  tryCatch(
    {
      if (!require(package, quietly = TRUE, character.only = TRUE)) {
        devtools::install_github(paste0('github::', package))
        message(paste0(package, ' successfully installed'))
      } else {
        message(paste0(package, ' already installed'))
      }
    },
    error = function(e){
      message(paste0('Error installing package ', package, '\nOriginal message:'))
      message(e)
    }#,
    # warning=function(w){
    #   message(paste0('Warning installing package ', lib, '\nOriginal message:'))
    #   message(w)
    # }
  )
}

renv::snapshot()
