#!/bin/bash -e
docker build -t m5c . # > /dev/null 2>&1
time docker run \
    -v ./ref:/asc25/ref \
    -v ./SRR23538290:/asc25/SRR23538290 \
    -v ./SRR23538291:/asc25/SRR23538291 \
    -v ./SRR23538292:/asc25/SRR23538292 \
    -it m5c > $(date +"%Y-%m-%d_%H-%M-%S").log