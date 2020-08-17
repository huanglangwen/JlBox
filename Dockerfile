FROM ubuntu:20.04
#MAINTAINER Langwen Huang (huanglangwen@outlook.com)
RUN mkdir -p /Code
RUN mkdir -p /Code/julia

WORKDIR /Code/julia
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y wget git && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.0-linux-x86_64.tar.gz && \
    tar xf julia-1.5.0-linux-x86_64.tar.gz && \
    rm julia-1.5.0-linux-x86_64.tar.gz && \
    echo "export PATH=/Code/julia/julia-1.5.0/bin:/root/.julia/conda/3/bin:$PATH" >> /root/.bashrc
ENV PYTHON=""
ENV PATH="/Code/julia/julia-1.5.0/bin:/root/.julia/conda/3/bin:${PATH}"
ADD . /Code/julia/JlBox
RUN julia --eval 'using Pkg;Pkg.develop(PackageSpec(path="/Code/julia/JlBox"));Pkg.build("JlBox")'

WORKDIR /root
RUN git clone git://github.com/JuliaEditorSupport/julia-vim.git
WORKDIR /root/julia-vim
RUN mkdir -p /root/.vim
RUN cp -R * /root/.vim
RUN echo "set nowrap\nset number\nset tabstop=4\nset shiftwidth=4\nset softtabstop=4\nset expandtab\nset smarttab" >> /root/.vimrc

WORKDIR /Code/julia/JlBox
RUN julia example/Install_package.jl