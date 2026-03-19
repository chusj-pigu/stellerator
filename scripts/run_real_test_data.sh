cargo run -- \
    --bam /home/geonic/Documents/GitHub/classy/test_data/200193_SISJ1435_T_asWGS_Sar_200245_1T_hg38.indexed.bam \
    --annotation /home/geonic/Downloads/sj-panel/hg38.ncbiRefSeq.gtf \
    --gene EWSR1 \
    --partner-gene FLI1 \
    --output-tsv test_data/output/stellerator.tsv \
    --output-fasta test_data/output/stellerator.fasta.gz \
    --threads 4 \
    --verbose
