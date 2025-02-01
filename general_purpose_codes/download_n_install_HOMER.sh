#!/usr/bin/env bash

#Download HOMER

wget http://homer.ucsd.edu/homer/configureHomer.pl 

#Install HOMER

./configureHomer.pl -install 

export PATH=$PATH:/data/collaborators/aak/NF1_sequencing/.//bin/ 
