#rm lists/jet30.imagelist
#rm lists/mb.imagelist
#rm lists/jet10.imagelist
#for rn in `ls  lists/dst_calo_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    ls /sphenix/tg/tg01/jets/jocl/evt/$rn/*alltog*$rn* > lists/$rn".imagelist"
#done
#
#for i in {0..9999}; do
#    echo $i
#    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*20250205_jet10* >> lists/jet10.imagelist
#done

#for i in {0..9999}; do
#    echo $i
#    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*20250205_jet30* >> lists/jet30.imagelist
#done
#exit
for i in {0..99999}; do
    echo $i
    ls /sphenix/tg/tg01/jets/jocl/evt/$i/*20250205_mb* >> lists/mb.imagelist
done
