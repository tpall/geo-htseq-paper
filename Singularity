BootStrap: shub
From: tpall/singularity-r:latest

%labels
  Maintainer tpall

%help
  This will run R packages to fit stan models

%post
  ## Download and install tidyverse & other packages
  apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    libssh2-1-dev \
    libudunits2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgdal-dev \
    libgsl-dev \
    libnode-dev \
    libnlopt-dev \
    libfontconfig1-dev \
    libmagick++-dev \
    cargo \
    jags \
    libharfbuzz-dev \
    libfribidi-dev \
  && Rscript -e 'install.packages(c("dplyr","readr","purrr","stringr","tidyr","lubridate","here","ggplot2","gt","extrafont","cowplot","patchwork","viridis","brms","rstan","tidybayes","magick","glue","BH","Rcpp", "RcppArmadillo"),type="source", dependencies=c("Depends", "Imports", "LinkingTo"))'
