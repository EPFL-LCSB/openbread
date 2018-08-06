#!/bin/bash

cd /src

git clone https://github.com/weilandtd/openfpm_pdata.git

cd openfpm_pdata

 ./install -s -i "/home/user/" \
           -c "--prefix=/usr/local/  CXXFLAGS=-fPIC  CFLAGS=-fPIC"


make install
