# README  

The runtime of each step in the workflow, along with the total execution time, is recorded in log files as follows:  

- `result_baseline.log`    # Baseline results  
- `result_baseline_ht.log` # Baseline results (with hyperthreading)  
- `result_final.log`       # Final optimized results  
- `result_final_ht.log`    # Final optimized results (with hyperthreading)  

Our best results are recorded in `result_final_ht.log`, while the other logs serve as comparisons.  

Our optimization of the workflow mainly focuses on process parallelism implemented in the shell and performance optimization for the hisat-3n project. To learn more, you can read the source code of run.sh and the README in the hisat-3n folder. Git commits are also a great way to explore, as our commits are well-organized and highly readable.

To quickly reproduce the results, use the `run_docker.sh` script, which will automatically build the image and run the container.  

If you want to install it locally, you can also refer to the configuration in the Dockerfile for the compilation and installation process.

```
.
├── bin
│   ├── filter_sites.py
│   ├── group_pileup.py
│   ├── join_pileup.py
│   └── select_sites.py
├── Dockerfile
├── environment.yml
├── fetch.sh               # Download data
├── hisat-3n/
├── htslib/
├── index.sh               # Indexing script
├── README.md
├── ref
│   ├── genome/
│   └── index/
├── result_baseline_ht.log # Baseline results (with hyperthreading)
├── result_baseline.log    # Baseline results
├── result_final_ht.log    # Final optimized results (with hyperthreading)
├── result_final.log       # Final optimized results
├── run_docker.sh          # Build and run container
├── run.sh                 # Workflow description file 
├── samtools/
├── SRR23538290
│   ├── SRR23538290.fastq
│   └── ...
├── SRR23538290.filtered.tsv
├── SRR23538291
│   └── ...
├── SRR23538291.filtered.tsv
├── SRR23538292
│   └── ...
├── SRR23538292.filtered.tsv
├── UMICollapse/
```