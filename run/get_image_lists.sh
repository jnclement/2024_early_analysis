for rn in `ls *.list | awk -F'.' '{print $1}'`; do
    rn=$(expr $rn + 0)
    ls output/evt/*newdat*$rn* > $rn".imagelist"
done

ls output/evt/*newsim* > sim.imagelist

#for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*0000* | awk -F'-' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    ls output/evt/*$rn* > $rn"_nszs.imagelist"
#done
