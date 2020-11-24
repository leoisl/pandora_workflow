# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_analysis_pipeline:pandora_paper_tag1

# This container has pandora to be used by the pandora_analysis_pipeline
# to build: sudo docker build . -t leandroishilima/pandora_analysis_pipeline:pandora_paper_tag1
FROM ubuntu:bionic

ARG pandora_version="pandora_paper_tag1"

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# install basic dependencies
RUN apt update && \
    apt install -y software-properties-common=0.96.24.32.14 && \
    apt-add-repository universe && \
    apt update && \
    apt install -y git wget \
        build-essential=12.4ubuntu1 \
        zlib1g-dev=1:1.2.11.dfsg-0ubuntu2 \
        python3=3.6.7-1~18.04 \
        python3-pip=9.0.1-2.3~ubuntu1.18.04.2 \
        cmake=3.10.2-1ubuntu2.18.04.1


# ========================
# INSTALL clustal omega
# ========================
RUN apt install -y clustalo=1.2.4-1


# ========================
# INSTALL boost_1_62_0
# ========================
RUN wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz -O - | tar xzf - && \
    cd boost_1_62_0 && \
    ./bootstrap.sh [--prefix=/prefix/path] --with-libraries=system,filesystem,iostreams,log,thread,date_time && \
    ./b2 install

# ========================
# INSTALL pandora
# ========================
RUN git clone https://github.com/rmcolq/pandora.git --recursive && \
    cd pandora && \
    git checkout $pandora_version && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=ReleaseWithAsserts .. && \
    make -j4 && \
    ctest -VV
ENV PATH="/pandora/build:${PATH}"
