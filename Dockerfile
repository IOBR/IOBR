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
    BiocManager::install(ask = FALSE, update = FALSE)"

# Install IOBR dependencies that are commonly needed
RUN Rscript -e "BiocManager::install(c( \
    'GSVA', \
    'limma', \
    'DESeq2', \
    'ComplexHeatmap', \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'preprocessCore', \
    'sva', \
    'biomaRt', \
    'enrichplot', \
    'DOSE' \
    ), ask = FALSE, update = FALSE)"

# Install additional CRAN packages
RUN Rscript -e "install.packages(c( \
    'survminer', \
    'glmnet', \
    'ggpubr', \
    'ggsci', \
    'RColorBrewer', \
    'scales', \
    'patchwork', \
    'corrplot', \
    'factoextra', \
    'FactoMineR', \
    'pROC', \
    'timeROC', \
    'tidyHeatmap', \
    'Hmisc', \
    'msigdbr', \
    'MASS', \
    'e1071', \
    'psych', \
    'reshape2', \
    'gridExtra', \
    'WGCNA' \
    ), repos = 'https://cloud.r-project.org')"

# Install IOBR from GitHub
RUN Rscript -e "BiocManager::install('IOBR/IOBR', ask = FALSE, update = FALSE, dependencies = TRUE)"

# Set working directory
WORKDIR /home/rstudio

# Set default command
CMD ["/init"]
