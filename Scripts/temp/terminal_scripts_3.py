mkdir -p project/task3/{lists,mhc_tables,logs}

./scripts/filter_to_mhc.py | tee project/task3/logs/task3_step1_filter.log

wc -l project/task3/mhc_tables/MOT36308.mhc_genes.tsv
head project/task3/mhc_tables/MOT36308.mhc_genes.tsv

