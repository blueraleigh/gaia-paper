for j in {1..25}; do
    Rscript --vanilla ooa.R $j north
    Rscript --vanilla ooa.R $j south
    Rscript --vanilla ooa.R $j both
done
