CWD=$(pwd)

for j in {1..25}; do
    slim -d WD=\"$CWD\" -d REP=$j -d MAP=\"north\" run.slim
    slim -d WD=\"$CWD\" -d REP=$j -d MAP=\"south\" run.slim
    slim -d WD=\"$CWD\" -d REP=$j -d MAP=\"both\" run.slim
done
