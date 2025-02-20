FROM intel/oneapi-hpckit:2025.0.2-0-devel-ubuntu24.04
WORKDIR /asc25
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /opt/conda \
    && rm -rf /tmp/miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN cat > ~/.condarc <<EOF
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
EOF

RUN conda install -y python=3 \
    && conda update conda \
    && conda clean --all --yes \
    && conda init

RUN apt-get -qq update \ 
    && apt-get -qq -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev bc \
    libcurl4-gnutls-dev libssl-dev libdeflate-dev libncurses-dev parallel cmake openjdk-11-jre numactl

COPY . .

RUN cd htslib && git checkout 1.21 && git submodule update --init --recursive \
    && make clean && autoreconf -i && ./configure CC=icx \
    && make -j$(nproc)
RUN cd samtools && git checkout 1.21 \
    && make clean && autoheader && autoconf -Wno-syntax && ./configure CC=icx \
    && make -j$(nproc)
RUN cd hisat-3n && git checkout optimized \
    && rm -rf build && cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -S . -B build \
    && cmake --build build -- -j$(nproc)

RUN pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
RUN conda env create -f environment.yml
CMD ["./run.sh"]
