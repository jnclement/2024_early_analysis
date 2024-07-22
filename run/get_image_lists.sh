rm lists/sim.imagelist

for rn in `ls lists/*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}'`; do
    rn=$(expr $rn + 0)
    ls /sphenix/tg/tg01/jets/jocl/evt/$rn/*newdat*$rn* > lists/$rn".imagelist"
done

for i in {0..20000}; do
    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*newsim* >> lists/sim.imagelist
done
