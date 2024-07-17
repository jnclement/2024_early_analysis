rm sim.imagelist

for rn in `ls *.list | awk -F'.' '{print $1}'`; do
    rn=$(expr $rn + 0)
    ls /sphenix/tg/tg01/jets/jocl/evt/$rn/*newdat*$rn* > $rn".imagelist"
done

for i in {0..10000}; do
    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*newsim* >> sim.imagelist
done
