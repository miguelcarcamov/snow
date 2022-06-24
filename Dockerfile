FROM ubuntu:20.04
RUN apt-get update -y && \
  DEBIAN_FRONTEND=noninteractive && \
  apt-get install -y tzdata && \
  apt-get install -y keyboard-configuration

RUN apt-get install -y build-essential && \
    apt-get install -y zlib1g-dev libncurses5-dev && \
    apt-get install -y libgdbm-dev libnss3-dev libssl-dev  && \
    apt-get install -y libreadline-dev libffi-dev wget && \
    apt-get install -y --no-install-recommends && \
    apt-get install -y python3-dev && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-wheel && \
    apt-get install -y python3-setuptools && \
    apt-get install -y libblas-dev && \
    apt-get install -y liblapack-dev && \
    apt-get install -y liblapacke-dev && \
    apt-get install -y git && \
    apt-get install -y ImageMagick* && \
    apt-get install -y xorg && \
    apt-get install -y compat-libgfortran-48 && \
    apt-get install -y libnsl && \
    apt-get install -y openmpi-devel && \
    apt-get install -y mpich-devel && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from selfcalframework base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/selfcalframework"
