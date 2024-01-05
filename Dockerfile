FROM ubuntu:jammy-20231211.1
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME) main restricted universe multiverse" > /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME) main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt $(. /etc/os-release && echo $VERSION_CODENAME)-security main restricted universe multiverse" >> /etc/apt/sources.list && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends ca-certificates && \
    apt-get install -y --no-install-recommends openssl && \
    update-ca-certificates && \
    echo 'APT::Get::Always-Include-Phased-Updates "true";' | tee -a /etc/apt/apt.conf.d/99-phased-updates && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends systemd libsystemd-dev && \
    apt-get install -y --no-install-recommends software-properties-common && \
    add-apt-repository main && \
    add-apt-repository universe && \
    add-apt-repository restricted && \
    add-apt-repository multiverse && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends tzdata && \
    apt-get install -y --no-install-recommends keyboard-configuration && \
    apt-get install -y --no-install-recommends build-essential && \
    apt-get install -y --no-install-recommends zlib1g-dev libncurses5-dev && \
    apt-get install -y --no-install-recommends libgdbm-dev libnss3-dev libssl-dev && \
    apt-get install -y --no-install-recommends libreadline-dev libffi-dev wget && \
    apt-get install -y --no-install-recommends --no-install-recommends && \
    apt-get install -y --no-install-recommends  python3-dev && \
    apt-get install -y --no-install-recommends python3-pip && \
    apt-get install -y --no-install-recommends python3-wheel && \
    apt-get install -y --no-install-recommends python3-setuptools && \
    apt-get install -y --no-install-recommends libblas-dev && \
    apt-get install -y --no-install-recommends liblapack-dev && \
    apt-get install -y --no-install-recommends liblapacke-dev && \
    apt-get install -y --no-install-recommends git && \
    apt-get install -y --no-install-recommends ImageMagick* && \
    apt-get install -y --no-install-recommends xorg && \
    apt-get install -y --no-install-recommends libgfortran4 && \
    apt-get install -y --no-install-recommends libopenmpi-dev && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN python3 --version
RUN pip3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install mpi4py --no-cache-dir
RUN echo "Hello from SNOW base image"
LABEL org.opencontainers.image.source="https://github.com/miguelcarcamov/snow"
LABEL org.opencontainers.image.description="Container image for SNOW"
LABEL org.opencontainers.image.licenses=GPL3
