#!/bin/bash -e
prefetch SRR23538290; prefetch SRR23538291; prefetch SRR23538292
fasterq-dump --threads 64 SRR23538290/SRR23538290.sra
fasterq-dump --threads 64 SRR23538291/SRR23538291.sra
fasterq-dump --threads 64 SRR23538292/SRR23538292.sra