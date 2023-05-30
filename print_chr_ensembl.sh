awk '{ if($1 !~ /^#/){print "chr"$0} else{print $0} }'  Mus_musculus.GRCm39t.106.chr.gtf > output.gtf
