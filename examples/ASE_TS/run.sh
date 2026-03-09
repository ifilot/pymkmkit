#!/bin/bash

# unpack archive
unzip ../../tests/data/CeO2_Pd4_COox.zip

# run pymkmkit
pymkmkit asevib2yaml OUTCAR -o ceo2_pd4_coox_ts.yaml