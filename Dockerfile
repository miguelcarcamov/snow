FROM ubuntu:20.04

RUN apt-get update -y && \
    apt-get install -y build-essential && \
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
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN pip3
RUN echo "Hello from selfcalframework base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/selfcalframework"
