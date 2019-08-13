FROM ubuntu:18.04
#MAINTAINER Langwen Huang (huanglangwen@outlook.com)
#forked from https://github.com/loftytopping/PyBox/blob/master/Dockerfile

RUN apt-get update 
RUN apt-get install -y build-essential \
    apt-utils \
    wget \
    vim \
    git \
    tmux \
    curl

RUN mkdir -p /Code
RUN mkdir -p /Code/julia

WORKDIR /Code/julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.1-linux-x86_64.tar.gz
RUN tar xf julia-1.1.1-linux-x86_64.tar.gz
RUN echo "export PATH=/Code/julia/julia-1.1.1/bin:/root/.julia/conda/3/bin:$PATH" >> /root/.bashrc
ENV PYTHON=""
ENV PATH="/Code/julia/julia-1.1.1/bin:/root/.julia/conda/3/bin:${PATH}"
#RUN source /root/.bashrc
RUN julia --eval 'using Pkg;Pkg.update()'
RUN julia --eval 'using Pkg;Pkg.add("Conda")'
RUN julia --eval 'using Conda;Conda.update()'
RUN julia --eval 'using Pkg;Pkg.add("PyCall")'
RUN julia --eval 'using Pkg;Pkg.build("PyCall")'
RUN julia --eval 'using Pkg;Pkg.dev("https://github.com/huanglangwen/JlBox");Pkg.build("JlBox")'

WORKDIR /root
RUN git clone git://github.com/JuliaEditorSupport/julia-vim.git
WORKDIR /root/julia-vim
RUN mkdir -p /root/.vim
RUN cp -R * /root/.vim
RUN echo "set nowrap\nset number\nset tabstop=4\nset shiftwidth=4\nset softtabstop=4\nset expandtab\nset smarttab" >> /root/.vimrc

WORKDIR /root/.julia/dev/JlBox
RUN git checkout master