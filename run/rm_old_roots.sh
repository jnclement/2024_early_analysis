for file in `ls output/temphists/*`; do 
    find $file -mtime +1 -exec rm {} \;
done
