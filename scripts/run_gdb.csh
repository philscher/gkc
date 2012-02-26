#!/bin/csh -f

echo "Running GDB on node `hostname`"
xterm -e gdb $*
exit 0

