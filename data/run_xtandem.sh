#!/bin/bash

# activate conda environment
source /opt/anaconda3/bin/activate xtandem_env

# go to project folder
cd ~/xtandem_test

# run X! Tandem
xtandem input.xml > log.txt 2>&1

# show output file info
ls -lh output.xml

echo "Done! Results saved in output.xml"
