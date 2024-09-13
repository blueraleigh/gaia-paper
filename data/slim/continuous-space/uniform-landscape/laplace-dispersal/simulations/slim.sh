SIGMA=(0.2 0.5 0.8 1.1 1.4 1.7 2.0)
CWD=$(pwd)

for i in ${!SIGMA[@]}; do
    for j in {0..9}; do
        slim -d S=${SIGMA[i]} -d WD=\"$CWD\" -d REP=$j run.slim
    done
done
