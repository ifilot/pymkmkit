#!/bin/bash

# unpack archive
unzip ../../tests/data/CeO2_Pd4_CO.zip

# run pymkmkit
pymkmkit asevib2yaml OUTCAR -o ceo2_pd4_co.yaml