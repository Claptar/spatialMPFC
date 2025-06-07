# 0. Base image: R 4.4.0 + full tidyverse already baked in
FROM rocker/tidyverse:4.4.0


# 1. Extra system libraries needed by Python scientific stack & Jupyter
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential git curl ca-certificates \
    libblas3 liblapack3 libatlas-base-dev libgl1 \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    nodejs && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# 2. Copy statically-linked uv binaries & install Python 3.11 system-wide
COPY --from=ghcr.io/astral-sh/uv:0.7.12 /uv /uvx /bin/
RUN uv python install 3.11


# 3. Create an isolated virtual environment in /env and prepend it to PATH
ENV VENV_PATH="/env"
ENV PATH="${VENV_PATH}/bin:${PATH}"
RUN uv venv "${VENV_PATH}"


# 4. Install all required Python packages with uv
RUN . "${VENV_PATH}/bin/activate" && \
    uv pip install scanpy squidpy tqdm jupyter gseapy decoupler \
    plotly seaborn ipywidgets notebook papermill


# 5. Add missing R & Bioconductor pieces, then register R with Jupyter
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com/')" && \
    R -e "BiocManager::install('edgeR', ask = FALSE, update = FALSE)" && \
    R -e "install.packages('IRkernel', repos='https://cran.rstudio.com/')" && \
    R -e "IRkernel::installspec(user = FALSE)"


# 6. Quality-of-life environment variables
ENV PYTHONUNBUFFERED=1 \
    TQDM_NOTEBOOK=1


# 7. Working directory, port & default command
WORKDIR /workspace
EXPOSE 8888


# 8. Copy Dockerfile to the container
COPY Dockerfile /docker/
RUN chmod -R 755 /docker

# 9. Set the default command to run Jupyter Lab
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root"]