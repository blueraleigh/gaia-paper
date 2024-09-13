for j in {1..20}; do
    Rscript --vanilla performance.R $j ''
done

for j in {1..20}; do
    Rscript --vanilla performance.R $j \-pareto
done
