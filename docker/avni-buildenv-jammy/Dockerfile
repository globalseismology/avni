FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    g++ \
    gfortran \
    git \
    libopenmpi-dev \
    make \
    nano \
    wget \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-funcsigs \
    python3-nbsphinx \
    jupyter \
    && apt-get clean && rm -r /var/lib/apt/lists/*

COPY requirements.txt requirements.txt

RUN pip3 install -r requirements.txt

# Export compilers
ENV F90 gfortran
