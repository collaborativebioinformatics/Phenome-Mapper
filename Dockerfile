FROM quay.io/vgteam/vg:v1.68.0

SHELL ["/bin/bash", "-lc"]
ARG DEBIAN_FRONTEND=noninteractive

# OS deps for Miniforge bootstrap
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates bzip2 \
 && rm -rf /var/lib/apt/lists/*

# Install Miniforge (arm64 and x86_64 builds available — adjust arch if needed)
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN curl -L -o Miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
 && bash Miniforge.sh -b -p $CONDA_DIR \
 && rm Miniforge.sh \
 && conda clean -afy


# install based on a recipe
# COPY conda_spec.txt /tmp/conda_spec.txt
RUN conda config --system --add channels conda-forge \
 && conda config --system --add channels bioconda \
 && conda config --system --set channel_priority flexible \
 && conda create --name phenomemapper spades quast minigraph \
	bcftools samtools pandas numpy seaborn \
	r-base r-data.table r-tidyverse \
	plink2 regenie \
 && conda clean -afy

# Sanity check (optional)
# RUN conda run -n phenomemapper spades.py --version && vg version

ENTRYPOINT ["conda", "run", "--live-stream", "-n", "phenomemapper"]
CMD ["bash"]
