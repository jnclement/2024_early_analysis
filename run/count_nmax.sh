rm listoflist.txt
for list in *.list; do
    cat $list | head -n 350 >> listoflist.txt
done
wc -l listoflist.txt
