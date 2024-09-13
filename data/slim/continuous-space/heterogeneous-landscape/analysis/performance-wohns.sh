for j in {1..20}; do
    python3 performance-wohns.py $j trees
done

for j in {1..20}; do
    python3 performance-wohns.py $j trees-pareto
done
