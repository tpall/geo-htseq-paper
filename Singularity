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
    git \
    curl \
  && Rscript -e 'install.packages(c("dplyr","readr","purrr","stringr","tidyr","lubridate","here","ggplot2","gt","extrafont","cowplot","patchwork","viridis","magick","glue","BH","Rcpp","RcppArmadillo","rmarkdown","bookdown","parallel","V8","webshot"),type="source", dependencies=c("Depends", "Imports", "LinkingTo"))' \
  && Rscript -e 'remotes::install_github("stan-dev/rstan",ref="develop",subdir="rstan/rstan",build_opts="")' \
  && Rscript -e 'install.packages(c("brms","tidybayes"))'

## C++ toolchain configuration
mkdir -p $HOME/.R \
  && printf "\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC\nCXX14=g++" > $HOME/.R/Makevars

## Clean up from R source install
  cd / \
  && apt-get autoremove -y \
  && apt-get autoclean -y \
  && rm -rf /var/lib/apt/lists/*

%test
    Rscript -e 'library(rstan); example(stan_model, package = "rstan", run.dontrun = TRUE)'

