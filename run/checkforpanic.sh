for rn in `cat runlist_forcheck.txt`; do
    root -l -b -q -e '.x check_scaledowns.C('${rn}')'
done
