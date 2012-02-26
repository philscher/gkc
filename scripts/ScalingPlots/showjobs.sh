qstat | egrep  -o [0-9]+.nqs | sed '{:q;N;s/\n/ /g;t q}'
