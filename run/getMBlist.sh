#!/bin/bash

#G4HITS_LIST="g4hits.list"
DST_CALO_CLUSTER_LIST="dst_calo_cluster.list"
DST_TRUTH_JET_LIST="dst_truth_jet.list"
#DST_GLOBAL_LIST="dst_global.list"
DST_MBD_EPD_LIST="dst_mbd_epd.list"

# Run SQL queries and redirect results to files
#psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 21 AND dsttype = 'G4Hits' AND filename LIKE '%Detroit%' ORDER BY segment;" > "$G4HITS_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 21 AND dsttype = 'DST_CALO_CLUSTER' AND filename LIKE '%Detroit%' ORDER BY segment;" > "$DST_CALO_CLUSTER_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 21 AND dsttype = 'DST_TRUTH_JET' AND filename LIKE '%Detroit%' ORDER BY segment;" > "$DST_TRUTH_JET_LIST"
#psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 21 AND dsttype = 'DST_GLOBAL' AND filename LIKE '%Detroit%' ORDER BY segment;" > "$DST_GLOBAL_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 21 AND dsttype = 'DST_MBD_EPD' AND filename LIKE '%Detroit%' ORDER BY segment;" > "$DST_MBD_EPD_LIST"
