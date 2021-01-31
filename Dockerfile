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

RUN apt-get update && apt-get install -y git jupyter-client jupyter-notebook && apt-get clean
RUN pip3 install --upgrade tornado jupyter jupyter_client notebook jupyterlab nbconvert

RUN mkdir -p ${HOME}/JlBox
COPY . ${HOME}/JlBox
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}/JlBox
RUN julia --eval 'using Pkg;Pkg.develop(PackageSpec(path="."));Pkg.build("JlBox")'
RUN julia example/Install_package.jl
RUN julia --eval 'using Pkg;Pkg.add("IJulia")'
RUN julia --eval 'using Pkg;Pkg.activate(".");Pkg.instantiate();using JlBox;using Plots'
# docker build . -t jlbox
# docker run --rm -p 8888:8888 jlbox jupyter notebook --ip=0.0.0.0 --no-browser
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()