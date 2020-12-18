BootStrap: shub
From: tpall/singularity-tidyverse:latest

%labels
  Maintainer tpall

%help
  This will run R packages to fit stan models

%post
  ## Download and install tidyverse & other packages
  apt-get update -qq \
  && apt-get -y --no-install-recommends install \
    libnlopt-dev \
    libfontconfig1-dev \
    libmagick++-dev \
    cargo \
    jags \
    libharfbuzz-dev \
    libfribidi-dev \
    git \
    curl \
    wget \
  && install2.r --error \
    --deps TRUE \
    --skipinstalled \
    brms \
    tidybayes \
    rstan \
    gt \
    extrafont \
    viridis \
    magick \
    V8 \
    sparkline

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

