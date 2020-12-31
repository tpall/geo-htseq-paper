BootStrap: docker
From: rstatstartu/rstanverse:latest

%labels
  Maintainer tpall

%help
  This will run R packages to fit stan models

%post
  Rscript -e "webshot::install_phantomjs()"

  cd / \
  && apt-get autoremove -y \
  && apt-get autoclean -y \
  && rm -rf /var/lib/apt/lists/*
