for j in {1..15}; do
    Rscript --vanilla performance-linear.R $j ''
done

for j in {1..15}; do
    Rscript --vanilla performance-linear.R $j \-pareto
done
