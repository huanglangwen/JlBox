FROM julia:1.5.2

ARG NB_USER=jlboxuser
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

RUN apt-get update && apt-get install -y git libffi-dev && apt-get clean

#RUN mkdir -p ${HOME}/JlBox
#COPY . ${HOME}/JlBox
RUN git clone https://github.com/huanglangwen/JlBox.git ${HOME}/JlBox
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}/JlBox
RUN julia --eval 'using Pkg;Pkg.develop(PackageSpec(path="."));Pkg.build("JlBox")'
RUN julia example/Install_package.jl
RUN julia --eval 'using Pkg;Pkg.activate(".");Pkg.instantiate();using JlBox;using Plots'
RUN julia --eval 'using Pkg;Pkg.activate(".");Pkg.instantiate();using Conda;Conda.pip_interop(true);Conda.pip("install","jupyterlab");Pkg.activate();ENV["JUPYTER"]=joinpath(Conda.BINDIR,"jupyter");Pkg.add("IJulia")'
# docker build . -t jlbox
# docker run --rm -it jlbox /bin/bash
# docker run --rm -p 8888:8888 jlbox jupyter notebook --ip=0.0.0.0 --no-browser
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()
