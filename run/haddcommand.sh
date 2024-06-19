if [ $# -lt 1 ]; then
    echo "need sim or dat (0 or 1)"
    exit 1
fi

TYPE=sim    

if [ $1 -eq 1 ]; then
    TYPE=dat
    for rn in `ls  *.imagelist | awk -F'.' '{print $1}'`; do
	hadd "output/sumroot/summed_${rn}_${TYPE}_base.root" `ls output/root/*$rn*${TYPE}*`
    done
    ls output/sumroot/*dat* > sumdatlist.list
fi



hadd "summed_${TYPE}.root" `ls output/root/*${TYPE}*`
