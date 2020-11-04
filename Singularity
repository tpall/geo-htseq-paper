BootStrap: shub
From: tpall/singularity-r:4.0.3

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
  && Rscript -e 'install.packages("devtools", dependencies = c("Depends", "Imports", "LinkingTo"))' \
  && Rscript -e 'devtools::install_github("r-rust/gifski@v0.8.3")' \
  && install2.r --error \
    --deps TRUE \
    --skipinstalled \
    readr \
    dplyr \
    ggplot2 \
    stringr \
    lubridate \
    purrr \
    gt \
    extrafont \
    cowplot \
    patchwork \
    viridisLite \
    brms \
    tidybayes \
    here
