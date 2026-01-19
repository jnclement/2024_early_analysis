

cd lists
if [ "$1" = "mb" ] || [ "$1" = "all" ]; then
    cd ../mblist
    CreateFileList.pl -nop -run 28 -type 26 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet5" ] || [ "$1" = "all" ]; then
    cd ../jet5list
    CreateFileList.pl -nop -run 28 -type 36 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet10" ] || [ "$1" = "all" ]; then
    cd ../jet10list
    CreateFileList.pl -nop -run 28 -type 12 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet15" ] || [ "$1" = "all" ]; then
    cd ../jet15list
    CreateFileList.pl -nop -run 28 -type 33 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hit
fi
if [ "$1" = "jet20" ] || [ "$1" = "all" ]; then
    cd ../jet20list
    CreateFileList.pl -nop -run 28 -type 21 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet30" ] || [ "$1" = "all" ]; then
    cd ../jet30list
    CreateFileList.pl -nop -run 28 -type 11 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet50" ] || [ "$1" = "all" ]; then
    cd ../jet50list
    CreateFileList.pl -nop -run 28 -type 34 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
if [ "$1" = "jet70" ] || [ "$1" = "all" ]; then
    cd ../jet70list
    CreateFileList.pl -nop -run 28 -type 35 DST_GLOBAL DST_CALO_CLUSTER DST_TRUTH_JET DST_MBD_EPD G4Hits
fi
#26 mb
#36 5 GeV
#12 10 GeV
#33 15 GeV
#21 20 GeV
#11 30 GeV
#34 50 GeV
#35 70 GeV
