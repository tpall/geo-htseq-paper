BootStrap: docker
From: rstatstartu/rstanverse:v0.2

%labels
  Maintainer tpall

%help
  This will run R packages to fit stan models

%post
  cd / \
  && apt-get autoremove -y \
  && apt-get autoclean -y \
  && rm -rf /var/lib/apt/lists/*

%test
  Rscript -e 'library(rstan); example(stan_model, package = "rstan", run.dontrun = TRUE)'
