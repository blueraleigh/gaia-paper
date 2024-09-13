SIGMA=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)

for s in ${SIGMA[@]}; do
    for j in {1..10}; do
        python3 performance-wohns.py $s $j
    done
done
