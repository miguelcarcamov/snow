FROM nvidia/cuda:12.1.1-devel-ubuntu22.04
ENV BUILD_DIR "build"
ENV CONFIG "Release"
ENV PATH /usr/local/cuda/bin${PATH:+:${PATH}}
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Update
RUN apt-get update -y

# Install python
RUN apt-get install -y  python3-dev && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-wheel && \
    apt-get install -y python3-setuptools

# Install dependencies
RUN apt-get install -y ca-certificates && \
    apt-get install -y openssl && \
    apt-get install -y software-properties-common && \
    update-ca-certificates && \
    apt-get update -y && \
    add-apt-repository main && \
    add-apt-repository universe && \
    add-apt-repository restricted && \
    add-apt-repository multiverse && \
    apt-get update -y

RUN apt-get install -y tzdata && \
    apt-get install -y keyboard-configuration && \
    apt-get install -y build-essential && \
    apt-get install -y cmake && \
    apt-get install -y gfortran && \
    apt-get install -y g++ && \
    apt-get install -y libncurses5-dev && \
    apt-get install -y libreadline-dev && \
    apt-get install -y flex && \
    apt-get install -y bison  && \
    apt-get install -y libblas-dev  && \
    apt-get install -y libcfitsio-dev  && \
    apt-get install -y wcslib-dev  && \
    apt-get install -y libfftw3-dev && \
    apt-get install -y libhdf5-serial-dev && \
    apt-get install -y python3-numpy && \
    apt-get install -y libboost-all-dev && \
    apt-get install -y libgsl-dev && \
    apt-get install -y zlib1g-dev libncurses-dev && \
    apt-get install -y libgdbm-dev libnss3-dev libssl-dev && \
    apt-get install -y libreadline-dev libffi-dev wget && \
    apt-get install -y --no-install-recommends && \
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

# Install casacore
RUN git clone --single-branch --branch v3.5.0 https://github.com/casacore/casacore.git && \
cd casacore && \
mkdir build && \
cd build && \
cmake -DUSE_FFTW3=ON -DUSE_OPENMP=ON -DUSE_HDF5=ON -DBUILD_PYTHON3=ON -DUSE_THREADS=ON -DBUILD_PYTHON=OFF .. && \
make -j2 && \
make install

# Install CUDA samples
RUN cd /usr/local/cuda && \
git clone --single-branch --branch v11.8 https://github.com/NVIDIA/cuda-samples.git samples && \
cd samples && \
mv Common common && \
mv Samples samples && \
cd common && \
mkdir inc && \
mv *.h inc/

# Install GPUVmem
RUN echo "Installing GPUVMEM"
RUN git clone https://github.com/miguelcarcamov/gpuvmem.git && \
    cd gpuvmem && \
    cmake . -B $BUILD_DIR -DCMAKE_BUILD_TYPE=$CONFIG && \
    cd $BUILD_DIR && \
    cmake --build . --target install --verbose -j `nproc`

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
