#!/bin/bash
#SBATCH --job-name=lct_autoindex
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jad33864@uga.edu

set -euo pipefail          # stop at the first error instead of continuing silently

ml vg/1.65.0

cd ~/pangenome/vg-builds
rm -f lct_index.*          # clear stale index files (your known autoindex gotcha)

vg autoindex --workflow sr-giraffe \
  -r ~/pangenome/ref/lct_enhancer.fa \
  -v ~/pangenome/ref/lct_enhancer_shifted.vcf.gz \
  -p lct_index \
  -t $SLURM_CPUS_PER_TASK \
  -V 2
