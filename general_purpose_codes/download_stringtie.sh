#!/usr/bin/env bash

#Download stringtie for assembly

git clone https://github.com/gpertea/stringtie
cd stringtie
make -j4 release
