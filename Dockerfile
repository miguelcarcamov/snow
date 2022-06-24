FROM ubuntu:20.04
RUN echo "deb mirror://mirrors.ubuntu.com/mirrors.txt focal main restricted universe multiverse" > /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt focal-updates main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt focal-security main restricted universe multiverse" >> /etc/apt/sources.list && \
    apt-get update -Y

RUN DEBIAN_FRONTEND=noninteractive && \
  apt-get install -y tzdata && \
  apt-get install -y keyboard-configuration

RUN add-apt-repository main && \
  add-apt-repository universe && \
  add-apt-repository restricted && \
  add-apt-repository multiverse && \
  apt-get update -y

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
    apt-get install -y libgfortran4 && \
    apt-get install -y libnsl-dev && \
    apt-get install -y libopenmpi-dev && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from selfcalframework base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/selfcalframework"
