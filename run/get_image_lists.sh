rm lists/sim.imagelist

#for rn in `ls lists/47522.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    ls /sphenix/tg/tg01/jets/jocl/evt/$rn/*newdat*$rn* > lists/$rn".imagelist"
#done

for i in {0..19999}; do
    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*detroit* >> lists/sim.imagelist
done
