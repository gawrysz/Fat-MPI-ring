#!/bin/bash

if [ $# != 1 ] ; then
	echo "Usage: $0 file_from_fatring"
	exit 1
fi

# The awk script selects the samplest with maximum data volume for the given run

gnuplot << EOF
se log xy
se xla "#chunks"
se yla "time (s)"
se ke le

p [][1e-5:100] "< awk 'BEGIN {n=-1} /ring test/ {n=\$8; print } {if (\$2*\$3 == n) print}' $1" u 2:4 t "setup", "" u 2:5 t "win create", "" u 2:6 t "win fence 1", "" u 2:7 t "Put | Get", "" u 2:8 t "win fence 2", "" u 2:8 t "win free", "" u 2:10 t "check",x*1e-4 t "1e-4 * #chunks"
pause mouse close
EOF
