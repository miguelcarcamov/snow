FROM ubuntu:jammy
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get update -y && \
    apt-get install -y ca-certificates && \
    apt-get install -y openssl && \
    apt-get install -y software-properties-common && \
    update-ca-certificates && \
    apt-get update -y && \
    add-apt-repository main && \
    add-apt-repository universe && \
    add-apt-repository restricted && \
    add-apt-repository multiverse && \
    apt-get update -y && \
    apt-get install -y tzdata && \
    apt-get install -y keyboard-configuration && \
    apt-get install -y build-essential && \
    apt-get install -y zlib1g-dev libncurses-dev && \
    apt-get install -y libgdbm-dev libnss3-dev libssl-dev && \
    apt-get install -y libreadline-dev libffi-dev wget && \
    apt-get install -y --no-install-recommends && \
    apt-get install -y  python3-dev && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-wheel && \
    apt-get install -y python3-setuptools && \
    apt-get install -y libblas-dev && \
    apt-get install -y liblapack-dev && \
    apt-get install -y liblapacke-dev && \
    apt-get install -y git && \
    apt-get install -y ImageMagick* && \
    apt-get install -y xorg && \
    apt-get install -y libgfortran5 && \
    apt-get install -y libopenmpi-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /tmp/* /var/tmp/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
