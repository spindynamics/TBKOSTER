FROM ubuntu:latest
MAINTAINER Pascal Thibaudeau <pthibaud@users.noreply.github.com>

ENV DEBIAN_FRONTEND noninteractive
ENV TZ=Europe/Paris

RUN \
    apt-get update && \
    apt-get install -y locales tzdata && \
    ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

#RUN \
#    locale-gen fr_FR.UTF-8 && \
#    dpkg-reconfigure locales && \
#    update-locale LANG="fr_FR.UTF-8" LANGUAGE="fr_FR"
#ENV LC_ALL="fr_FR.UTF-8"

RUN \
  apt-get -y install \
          build-essential \
          gfortran \
          libblas-dev \
          liblapack-dev \
          cmake \
          vim \
          doxygen \
          graphviz \
          texlive \
          texlive-latex-recommended \
          texlive-latex-extra \
          gettext \ 
          tex4ht \
          biber \ 
          ghostscript \
          poppler-utils \
          git-core && \
  rm -rf /var/lib/apt/lists/*

ARG USERNAME=TBKOSTER
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the local user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

USER $USERNAME
WORKDIR /home/$USERNAME

# clone the repo and build the executables
RUN git clone https://github.com/spindynamics/TBKOSTER.git
WORKDIR /home/$USERNAME/TBKOSTER
RUN \
    mkdir build && \
    cmake -B build . && \
    cd build && \
    make && \
    make doc