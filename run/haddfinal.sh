ls output/sumroot/*dat* > sumdatlist.list

for TYPE in sim dat; do
    hadd "summed_${TYPE}.root" `ls output/root/*${TYPE}*fullfile*`
done
