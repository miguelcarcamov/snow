FROM ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive
ENV TZ=Etc/UTC
RUN echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME) main restricted universe multiverse" > /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME) main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME)-security main restricted universe multiverse" >> /etc/apt/sources.list && \
    apt-get update -y && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* && \
    apt-get install -y keyboard-configuration --no-install-recommends && \
    apt-get install -y software-properties-common --no-install-recommends

RUN add-apt-repository main && \
  add-apt-repository universe && \
  add-apt-repository restricted && \
  add-apt-repository multiverse && \
  apt-get update -y && \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* &&\
  apt-get install -y build-essential --no-install-recommends && \
  apt-get install -y zlib1g-dev libncurses5-dev --no-install-recommends && \
  apt-get install -y libgdbm-dev libnss3-dev libssl-dev  --no-install-recommends && \
  apt-get install -y libreadline-dev libffi-dev wget --no-install-recommends && \
  apt-get install -y --no-install-recommends --no-install-recommends && \
  apt-get install -y python3-dev --no-install-recommends && \
  apt-get install -y python3-pip --no-install-recommends && \
  apt-get install -y python3-wheel --no-install-recommends && \
  apt-get install -y python3-setuptools --no-install-recommends && \
  apt-get install -y libblas-dev --no-install-recommends && \
  apt-get install -y liblapack-dev --no-install-recommends && \
  apt-get install -y liblapacke-dev --no-install-recommends && \
  apt-get install -y git --no-install-recommends && \
  apt-get install -y ImageMagick* --no-install-recommends && \
  apt-get install -y xorg --no-install-recommends && \
  apt-get install -y libgfortran4 --no-install-recommends && \
  apt-get install -y libopenmpi-dev --no-install-recommends && \
  rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
