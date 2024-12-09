rm lists/sumdatlist.list
for rn in `ls lists/*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}'`; do
    echo $rn
    ls output/sumroot2/*$rn*dat* >> lists/sumdatlist.list
done
