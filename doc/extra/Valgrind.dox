cat gkc_log | ../../scripts/parse_valgrind_supp.sh > gkc.supp
valgrind   --log-file=gkc_log --gen-suppressions=all ../../src/gkc -c Test.gkc -f 
