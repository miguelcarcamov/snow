FROM ghcr.io/miguelcarcamov/gpuvmem:latest

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
