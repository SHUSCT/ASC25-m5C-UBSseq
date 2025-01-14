FROM intel/oneapi-hpckit:2025.0.2-0-devel-ubuntu24.04
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

WORKDIR /asc25

COPY . .

RUN apt-get -qq update \ 
    && apt-get -qq -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev libdeflate-dev libncurses-dev parallel cmake openjdk-11-jre numactl
RUN cd htslib && git checkout 1.21 && git submodule update --init --recursive \
    && make clean && autoreconf -i && ./configure CC=icx \
    && make -j$(nproc)
RUN cd samtools && git checkout 1.21 \
    && make clean && autoheader && autoconf -Wno-syntax && ./configure CC=icx \
    && make -j$(nproc)
RUN cd hisat-3n \
    && rm -rf build && cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -S . -B build \
    && cmake --build build -- -j$(nproc)

RUN conda env create -f environment.yml
SHELL ["conda", "run", "--no-capture-output", "-n", "m5C", "/bin/bash", "-c"]

CMD ["./run.sh"]
