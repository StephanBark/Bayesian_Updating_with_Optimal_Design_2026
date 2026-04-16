# =============================================================================
# Dockerfile for Bayesian Updating with Optimal Design (R/Stan/Gurobi project)
# Shared base plus dedicated runtime targets for VS Code and RStudio Server
# =============================================================================

FROM rocker/verse:4.4.2 AS workspace-base

LABEL maintainer="Stephan Bark (stephan-bark@web.de)"
LABEL author="Stephan Bark (stephan-bark@web.de)"
LABEL description="Bayesian Updating with Optimal Design using R/Stan/Gurobi in a reusable workspace image"
LABEL version="1.0"
LABEL source="https://github.com/StephanBark/Bayesian_Updating_with_Optimal_Design_2026"
LABEL created="2024-06-01"

# -- Build argument: Gurobi version (override version with --build-arg) -----
ARG GUROBI_VERSION=13.0.1

# -- System dependencies for R packages --------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libglpk-dev \
    python3-pip \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# -- Configure C++ toolchain for Stan -----------------------------------------
RUN mkdir -p /home/rstudio/.R && \
    echo "CXX14FLAGS = -O3 -march=native -mtune=native -fPIC" \
         > /home/rstudio/.R/Makevars && \
    echo "CXX14 = g++"  >> /home/rstudio/.R/Makevars && \
    echo "CXX17FLAGS = -O3 -march=native -mtune=native -fPIC" \
         >> /home/rstudio/.R/Makevars && \
    echo "CXX17 = g++"  >> /home/rstudio/.R/Makevars && \
    chown -R rstudio:rstudio /home/rstudio/.R

# -- Install CRAN R packages --------------------------------------------------
RUN install2.r --error --skipinstalled --ncpus -1 \
    data.table \
    readxl \
    writexl \
    reshape2 \
    lubridate \
    fitdistrplus \
    matrixcalc \
    viridis \
    showtext \
    gridExtra \
    patchwork \
    GGally \
    MCMCpack \
    coda \
    ggmcmc \
    SpATS \
    OptimalDesign \
    slam \
    htmltools \
    base64enc \
    languageserver \
    && rm -rf /tmp/downloaded_packages

# -- Install rstan (Stan interface) -------------------------------------------
RUN install2.r --error --skipinstalled --ncpus -1 \
    StanHeaders \
    rstan \
    && rm -rf /tmp/downloaded_packages

# -- Pre-compile a trivial Stan model to cache the toolchain -------------------
RUN Rscript -e ' \
    library(rstan); \
    rstan_options(auto_write = TRUE); \
    options(mc.cores = parallel::detectCores()); \
    stancode <- "parameters { real y; } model { y ~ normal(0,1); }"; \
    mod <- stan_model(model_code = stancode); \
    message("Stan toolchain verified OK") \
'

# -- Gurobi Optimizer (required by OptimalDesign for od_MISOCP) ---------------
# The installer archive must be placed in the project root before building.
# Download from: https://www.gurobi.com/downloads/gurobi-software/
#
# To build WITH Gurobi:
#   1. Place gurobi<VERSION>_linux64.tar.gz in the project root
#   2. docker compose up --build -d
#
# To build WITHOUT Gurobi:
#   Simply build without placing the archive.
#   The image still builds, but the od_MISOCP-based optimal-design workflows in
#   this repository will not run until the gurobi R package is installed.
# -----------------------------------------------------------------------------
# The Dockerfil[e] glob always matches, so COPY never fails even when no
# gurobi archive is present.  The RUN step below checks whether the tar.gz
# actually exists before attempting to install.
COPY Dockerfil[e] gurobi*_linux64.tar.gz /tmp/
RUN if ls /tmp/gurobi*_linux64.tar.gz 1>/dev/null 2>&1; then \
      if [ ! -f "/tmp/gurobi${GUROBI_VERSION}_linux64.tar.gz" ]; then \
        FOUND=$(ls /tmp/gurobi*_linux64.tar.gz | head -1) && \
        echo "ERROR: Found Gurobi installer '${FOUND}' but GUROBI_VERSION is set to ${GUROBI_VERSION}." && \
        echo "       Place gurobi${GUROBI_VERSION}_linux64.tar.gz in the project root, or pass" && \
        echo "       --build-arg GUROBI_VERSION=<version> matching the archive you provided." && \
        exit 1; \
      fi && \
      mv "/tmp/gurobi${GUROBI_VERSION}_linux64.tar.gz" /tmp/gurobi_installer.tar.gz; \
    fi && \
    if [ -f /tmp/gurobi_installer.tar.gz ]; then \
      echo ">>> Installing Gurobi ${GUROBI_VERSION} ..." && \
      tar -xzf /tmp/gurobi_installer.tar.gz -C /opt && \
      GUROBI_SHORT=$(printf '%s' "${GUROBI_VERSION}" | tr -d '.') && \
      GUROBI_MAJOR_MINOR=$(printf '%s' "${GUROBI_VERSION}" | sed 's/\.[^.]*$//') && \
      GUROBI_HOME="/opt/gurobi${GUROBI_SHORT}/linux64" && \
      echo "export GUROBI_HOME=${GUROBI_HOME}" >> /etc/profile.d/gurobi.sh && \
      echo "export PATH=\${GUROBI_HOME}/bin:\${PATH}" >> /etc/profile.d/gurobi.sh && \
      echo "export LD_LIBRARY_PATH=\${GUROBI_HOME}/lib:\${LD_LIBRARY_PATH}" >> /etc/profile.d/gurobi.sh && \
      chmod +x /etc/profile.d/gurobi.sh && \
      echo "GUROBI_HOME=${GUROBI_HOME}" >> /usr/local/lib/R/etc/Renviron && \
      echo "LD_LIBRARY_PATH=${GUROBI_HOME}/lib" >> /usr/local/lib/R/etc/Renviron && \
      ln -sf "${GUROBI_HOME}/lib/libgurobi$(printf '%s' "${GUROBI_MAJOR_MINOR}" | tr -d '.').so" /usr/local/lib/ && \
      ldconfig && \
      R CMD INSTALL "${GUROBI_HOME}"/R/gurobi_*.tar.gz && \
      rm /tmp/gurobi_installer.tar.gz && \
      echo ">>> Gurobi ${GUROBI_VERSION} installed successfully"; \
    else \
      echo ">>> No Gurobi installer found – skipping Gurobi installation; od_MISOCP-based optimal-design workflows will remain unavailable"; \
    fi

# -- Asreml requires a commercial licence -------------------------------------
# With available licence, uncomment and adjust the line below:
# RUN R -e "install.packages('asreml', repos = 'https://vsni.r-universe.dev')"

# -- Radian as modern R console specially when working in VS Code --------------
RUN pip3 install --no-cache-dir --break-system-packages radian

# -- Project files -------------------------------------------------------------
WORKDIR /home/rstudio/project

COPY . /home/rstudio/project/

# Copy the CMU Serif font into the system fonts directory for figures ---------
RUN if [ -f cmunrm.ttf ]; then \
      mkdir -p /usr/share/fonts/truetype/cmu && \
      cp cmunrm.ttf /usr/share/fonts/truetype/cmu/ && \
      fc-cache -fv; \
    fi

RUN chown -R rstudio:rstudio /home/rstudio/project

FROM workspace-base AS devcontainer

LABEL description="Bayesian Updating with Optimal Design for VS Code Dev Containers"

CMD ["sleep", "infinity"]

FROM workspace-base AS rstudio

LABEL description="Bayesian Updating with Optimal Design in a Docker container with RStudio Server"

# -- Expose RStudio Server (default port) -------------------------------------
EXPOSE 8787

# -- Default: start RStudio Server --------------------------------------------
CMD ["/init"]
