FROM ubuntu:20.04

RUN apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install -y zlib1g-dev libncurses5-dev && \
    apt-get install -y libgdbm-dev libnss3-dev libssl-dev  && \
    apt-get install -y libreadline-dev libffi-dev wget && \
    apt-get install -y --no-install-recommends && \
    apt-get install -y libblas-dev && \
    apt-get install -y liblapack-dev && \
    apt-get install -y liblapacke-dev && \
    apt-get install -y git && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python --version
RUN pip --version
RUN python -m pip install --upgrade pip
RUN pip
RUN echo "Hello from selfcalframework base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/selfcalframework"
