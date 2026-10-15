FROM nvidia/cuda:12.1.1-devel-ubuntu22.04
# GPUVmem, and the casacore and CUDA samples builds it needed, were removed: Pyralysis replaces it
# as snow's GPU imager (#24). casatools and python-casacore ship their own casacore.
ENV PATH /usr/local/cuda/bin${PATH:+:${PATH}}
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Update and install all system dependencies in a single layer to reduce image size
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        python3-dev python3-pip python3-wheel python3-setuptools \
        ca-certificates openssl software-properties-common \
        tzdata keyboard-configuration \
        build-essential cmake gfortran g++ \
        libncurses5-dev libreadline-dev flex bison \
        libblas-dev liblapack-dev liblapacke-dev \
        libcfitsio-dev wcslib-dev libfftw3-dev libhdf5-serial-dev \
        python3-numpy libboost-all-dev libgsl-dev \
        zlib1g-dev libncurses-dev libgdbm-dev libnss3-dev libssl-dev \
        libreadline-dev libffi-dev wget git \
        libgfortran5 libopenmpi-dev \
    && update-ca-certificates \
    && add-apt-repository -y main universe restricted multiverse \
    && apt-get update -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* /tmp/* /var/tmp/*

# Upgrade pip and install Python dependencies (no cache to save space)
RUN python3 --version && \
    pip3 --version && \
    python3 -m pip install --upgrade --no-cache-dir pip && \
    python3 -m pip install --no-cache-dir mpi4py && \
    python3 -m pip cache purge && \
    echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
