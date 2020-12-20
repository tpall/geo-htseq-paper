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
    libharfbuzz-dev \
    libfribidi-dev \
    libavfilter-dev \
    libpoppler-cpp-dev \
    libtesseract-dev \
    libleptonica-dev \
    tesseract-ocr-eng \
    cargo \
    jags \
    g++ \
    git \
    curl \
    wget

# Install cmdstan
VERSION=${VERSION:-2.25.0}
curl -O https://github.com/stan-dev/cmdstan/releases/download/v${VERSION}/cmdstan-${VERSION}.tar.gz \
    && tar -xzvf cmdstan-${VERSION}.tar.gz \
    && ln -s cmdstan-${VERSION} cmdstan \
    && cd cmdstan \
    && make build

# We need to update repos to install cmdstanr
CRAN=$(Rscript -e 'cat(getOption("repos"))')

install2.r --deps TRUE \
    --skipinstalled \
    --repos $CRAN https://mc-stan.org/r-packages/ \
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
