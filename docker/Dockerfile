FROM python:3.5
# Install missing deps
USER root

RUN apt-get update && apt-get install -y  \
        libbz2-dev          \
        libxml2-dev         \
        libxslt1-dev        \
        libpng-dev          \
        zlib1g-dev          \
        libgmp3-dev         \
        libc6-dev-i386      \
        gfortran            \
        git                 \
        emacs			    \
        python-dev          \
        swig                \
        cmake               \
        autoconf            \
        autotools-dev       \
        libtool             \
        wget                \
        htop                \
    && rm -rf /var/lib/apt/lists/*

ENV USER user
ENV HOME /home/user


RUN useradd -ms "/bin/bash" "$USER"
WORKDIR $HOME

COPY requirements.txt .


USER root

# Make an alias
RUN ln -s /usr/bin/python3 python
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

COPY install_openfpm.sh .
RUN /bin/bash ./install_openfpm.sh

RUN cp /usr/local/openfpm_pdata/include/example.mk $HOME/example.mk
RUN chmod -R o+rwx $HOME/example.mk

# Make the workspace folder that will link the sources
RUN mkdir /pymes

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"

# Fix Matplotlib
COPY matplotlibrc /usr/local/lib/python3.5/site-packages/matplotlib/mpl-data/matplotlibrc

#Install patchelf
RUN  cd /tmp && wget http://nixos.org/releases/patchelf/patchelf-0.8/patchelf-0.8.tar.gz && tar xf patchelf-0.8.tar.gz && patchelf-0.8/configure && make install && rm -rf patchelf-0.8 && rm -f patchelf-0.8.tar.gz

#Finalizing installation
USER $USER
RUN mkdir $HOME/work
WORKDIR  $HOME/work


WORKDIR $HOME/work
USER user

ENV PYTHON_INCLUDE='/usr/local/include/python3.5m'


ENTRYPOINT ["/bin/bash", "-c",  "source /home/user/openfpm_vars \
                                && pip install --user -e /openbread  && $0 $*"]

CMD /bin/bash
