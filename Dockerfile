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

RUN apt-get update && apt-get install -y git jupyter-client && apt-get clean 

RUN mkdir -p ${HOME}/JlBox
COPY . ${HOME}/JlBox
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}/JlBox
RUN julia --eval 'using Pkg;Pkg.develop(PackageSpec(path="${HOME}/JlBox"));Pkg.build("JlBox")'
RUN julia example/Install_package.jl
RUN julia --eval 'using Pkg;Pkg.add("IJulia")'