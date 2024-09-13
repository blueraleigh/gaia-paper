CWD=$(pwd)

for j in {1..20}; do
    slim -d WD=\"$CWD\" -d REP=$j run.slim
    slim -d WD=\"$CWD\" -d REP=$j run-pareto.slim
done
