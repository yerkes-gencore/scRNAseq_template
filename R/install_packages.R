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

for (package in cran_packages){
  tryCatch(
    {
      if (!require(package, quietly = FALSE)){
        renv::install.packages(package)
      }
    },
    error=function(e){
      message(paste0('Error installing package ', lib, '\nOriginal message:'))
      message(e)
    },
    warning=function(w){
      message(paste0('Warning installing package ', lib, '\nOriginal message:'))
      message(w)
    }
  )
}

for (package in bioc_packages){
  tryCatch(
    {
      if (!require(package, quietly = FALSE)){
        renv::install(paste0('bioc::', package))
      }
    },
    error=function(e){
      message(paste0('Error installing package ', lib, '\nOriginal message:'))
      message(e)
    },
    warning=function(w){
      message(paste0('Warning installing package ', lib, '\nOriginal message:'))
      message(w)
    }
  )
}

for (package in github_packages){
  tryCatch(
    {
      if (!require(package, quietly = FALSE)){
        renv::install(paste0('github::', package))
      }
    },
    error=function(e){
      message(paste0('Error installing package ', lib, '\nOriginal message:'))
      message(e)
    },
    warning=function(w){
      message(paste0('Warning installing package ', lib, '\nOriginal message:'))
      message(w)
    }
  )
}

renv::snapshot()
