FROM ubuntu:18.04
#MAINTAINER Langwen Huang (huanglangwen@outlook.com)
#forked from https://github.com/loftytopping/PyBox/blob/master/Dockerfile

RUN apt-get update 
RUN apt-get install -y build-essential \
    apt-utils \
    wget \
    libgl1-mesa-glx \
    gfortran \
    linux-headers-generic \
    cmake \
    vim \
    unzip \
    pkgconf \
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1 \
    libxrender1 \
    libxml2 \
    xauth \
    git \
    tmux

RUN mkdir -p /Code
RUN mkdir -p /Code/julia
RUN mkdir -p /Code/Git_repos

WORKDIR /Code/julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz
RUN tar xf julia-0.6.4-linux-x86_64.tar.gz
RUN mv julia-9d11f62bcb julia-0.6.4
RUN echo "export PATH=/Code/julia/julia-0.6.4/bin:/root/.julia/v0.6/Conda/deps/usr/bin:$PATH" >> /root/.bashrc
ENV PYTHON=""
ENV PATH="/Code/julia/julia-0.6.4/bin:/root/.julia/v0.6/Conda/deps/usr/bin:${PATH}"
#RUN source /root/.bashrc
RUN julia --eval 'Pkg.update()'
RUN julia --eval 'Pkg.add("DifferentialEquations")'
RUN julia --eval 'Pkg.add("StaticArrays")'
RUN julia --eval 'Pkg.add("DataFrames")'
RUN julia --eval 'Pkg.add("CSV")'
RUN julia --eval 'Pkg.add("Conda");using Conda;Conda.update()'
RUN julia --eval 'Pkg.add("PyCall")'
RUN julia --eval 'Pkg.build("PyCall")'
RUN julia --eval 'Pkg.add("LightXML")'
RUN conda config --append channels conda-forge
RUN conda install -c openbabel -y openbabel 
RUN conda install -y flask flask-wtf xlsxwriter

WORKDIR /Code/Git_repos
RUN git clone https://github.com/loftytopping/PyBox.git
RUN git clone https://github.com/loftytopping/UManSysProp_public.git
RUN git clone https://github.com/huanglangwen/JlBox.git

WORKDIR /Code/Git_repos/JlBox
RUN git checkout oldversion

WORKDIR /root
RUN git clone git://github.com/JuliaEditorSupport/julia-vim.git
WORKDIR /root/julia-vim
RUN mkdir -p /root/.vim
RUN cp -R * /root/.vim