FROM ubuntu:bionic

# to avoid interaction with tzdata during the installation
ARG DEBIAN_FRONTEND=noninteractive

# to avoid the UnicodeEncodeError when running phylosofs -h (default LANG is C)
ENV LANG C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    python3 \
    python3-pip \
    # PyCall.jl needs libpython3.x
    libpython3.6 \
    # setuptools is needed to install phylosofs
    python3-setuptools \
    # wheel is needed to install networkx
    python3-wheel \
    # ca-certificates is needed to use git clone without SSL CA cert problems
    ca-certificates \
    # PhyloSofS dependencies:
    graphviz \
    # make and cmake are needed to install HH-suite
    cmake \
    make \
    # install cmake compilers
    gcc \
    g++ \
    # to avoid xdd not found in HH-suite compilation
    xxd \
    # install MPI for HH-suite
    libopenmpi2 \
    libblacs-mpi-dev \
    # curl is needed to download MODELLER and Julia
    curl \
    # zlib in needed to use BioAlignments.jl
    zlib1g-dev \
    # qt5 is needed to use GR.jl
    qt5-default \
    # X11, OpenGL and Xvfb to use GR (gr-framework.org/tutorials/docker.html)
    libxt6 \
    libxrender1 \
    libxext6 \
    libgl1-mesa-glx \
    xvfb

# Output to PDF images (gr-framework.org/tutorials/docker.html)
ENV GKS_WSTYPE=pdf

WORKDIR /app

RUN git clone https://github.com/PhyloSofS-Team/PhyloSofS.git && \
    python3 -m pip install ./PhyloSofS

RUN git clone https://github.com/AntoineLabeeuw/hh-suite.git && \
   mkdir -p hh-suite/build && cd hh-suite/build && \
   cmake -DCMAKE_INSTALL_PREFIX=. .. && \
   make -j 4 && make install

RUN git clone https://github.com/soedinglab/pdbx.git && \
    mkdir -p pdbx/build && cd pdbx/build && \
    cmake -DUserInstallOption=ON ../ && \
    make install

ENV PATH="/app/hh-suite/build/bin:/app/hh-suite/build/scripts:${PATH}"

RUN curl -O https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.1-linux-x86_64.tar.gz && \
    tar xvzf julia-1.1.1-linux-x86_64.tar.gz && \
    rm julia-1.1.1-linux-x86_64.tar.gz && \
    /app/julia-1.1.1/bin/julia -e \
    'using Pkg; Pkg.activate("/usr/local/lib/python3.6/dist-packages/phylosofs/src/"); Pkg.instantiate(); pkg"precompile"'

ENV PATH="/app/julia-1.1.1/bin:${PATH}"

RUN curl -O https://salilab.org/modeller/9.21/modeller_9.21-1_amd64.deb && \
    dpkg -i modeller_9.21-1_amd64.deb && \
    rm modeller_9.21-1_amd64.deb

COPY docker_banner.sh docker_banner.sh

RUN cat /app/docker_banner.sh >> ~/.bashrc

WORKDIR /databases

WORKDIR /project

CMD ["/bin/bash"]
