FROM rocker/tidyverse:latest

LABEL org.opencontainers.image.source="https://github.com/IOBR/IOBR"
LABEL org.opencontainers.image.description="IOBR: Immuno-Oncology Biological Research"
LABEL org.opencontainers.image.authors="IOBR Team"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgdal-dev \
    libproj-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and set up Bioconductor
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    install.packages('pak', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(ask = FALSE, update = FALSE);"

# Install IOBR from GitHub
RUN Rscript -e "pak::pkg_install('IOBR/IOBR', dependencies = TRUE)"

# Set working directory
WORKDIR /home/rstudio

# Set default command
CMD ["/init"]
