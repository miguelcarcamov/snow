FROM nvidia/cuda:12.1.1-devel-ubuntu22.04
ENV BUILD_DIR "build"
ENV CONFIG "Release"
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

# Install casacore and clean up build artifacts to save space
RUN git clone --single-branch --branch v3.5.0 --depth 1 https://github.com/casacore/casacore.git && \
    cd casacore && \
    mkdir build && \
    cd build && \
    cmake -DUSE_FFTW3=ON -DUSE_OPENMP=ON -DUSE_HDF5=ON -DBUILD_PYTHON3=ON -DUSE_THREADS=ON -DBUILD_PYTHON=OFF .. && \
    make -j$(nproc) && \
    make install && \
    cd / && \
    rm -rf /casacore && \
    rm -rf /tmp/* /var/tmp/*

# Install CUDA samples (shallow clone to save space)
RUN cd /usr/local/cuda && \
    git clone --single-branch --branch v11.8 --depth 1 https://github.com/NVIDIA/cuda-samples.git samples && \
    cd samples && \
    mv Common common && \
    mv Samples samples && \
    cd common && \
    mkdir inc && \
    mv *.h inc/ && \
    cd / && \
    rm -rf /tmp/* /var/tmp/*

# Install GPUVmem and clean up build artifacts
RUN git clone --depth 1 https://github.com/miguelcarcamov/gpuvmem.git && \
    cd gpuvmem && \
    cmake . -B $BUILD_DIR -DCMAKE_BUILD_TYPE=$CONFIG && \
    cd $BUILD_DIR && \
    cmake --build . --target install --verbose -j$(nproc) && \
    cd / && \
    rm -rf /gpuvmem && \
    rm -rf /tmp/* /var/tmp/*

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
