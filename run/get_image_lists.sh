for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/temp24DSTs/*0000* | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    ls output/evt/*$rn* > $rn"_yszs.imagelist"
done

for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*0000* | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    ls output/evt/*$rn* $rn"_nszs.imagelist"
done
