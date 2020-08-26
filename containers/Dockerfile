# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_analysis_pipeline:pandora_paper_tag1

# This container has pandora and make_prg_dependencies to be used by the pandora_analysis_pipeline
# to build: sudo docker build . -t leandroishilima/pandora_analysis_pipeline:pandora_paper_tag1
FROM ubuntu:bionic

ARG pandora_version="pandora_paper_tag1"

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# install basic dependencies
RUN apt update && \
    apt install -y software-properties-common && \
    apt-add-repository universe && \
    apt update && \
    apt install -y git wget build-essential zlib1g-dev python3 python3-pip cmake


# ========================
# INSTALL clustal omega
# ========================
RUN apt install -y clustalo


# ========================
# INSTALL make_prg dependencies
# ========================
RUN pip3 install biopython==1.72 numpy==1.14.5 scikit-learn==0.20.3 snakemake && \
    cp /usr/bin/python3 /usr/bin/python


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
