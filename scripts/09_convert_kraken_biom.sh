#!/bin/bash
# scripts/09_kraken_biom_convert.sh

mamba run -n kraken-biom kraken-biom results/kraken2/*_report_bracken_species.tsv \
    -o results/bracken/combined.biom \
    --fmt json