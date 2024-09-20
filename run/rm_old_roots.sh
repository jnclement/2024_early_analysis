for file in `ls output/root/*`; do 
    find $file -mtime +2 -exec rm {} \;
done
