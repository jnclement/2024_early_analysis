hadd -j 4 -f /sphenix/user/jocl/projects/run2024_earlydata/run/output/chi2_haddhists/chi2_hadded_$1.root `cat /sphenix/user/jocl/projects/run2024_earlydata/run/chi2_hadd_condor_files.txt | head -n $(( $(( $1 + 1 )) * 100 )) | tail -n 100`