rm lists/sim.imagelist

#for rn in `ls  lists/dst_calo_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    ls /sphenix/tg/tg01/jets/jocl/evt/$rn/*alltog*$rn* > lists/$rn".imagelist"
#done

for i in {0..9999}; do
    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*goodcuts* >> lists/sim.imagelist
done
