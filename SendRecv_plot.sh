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

a=1
b=1e-10
c=1e-10
f(x) = a + b*x**2 + c*x**3
fit f(x) "< awk 'BEGIN {n=-1} /ring test/ {n=\$8; print } {if (\$2*\$3 == n) print}' $1" u 2:6:(column(2)) yerror via a,b,c

p [][1e-5:100] "" u 2:4 t "setup", "" u 2:5 t "SendRecv", "" u 2:6 t "Waitall", "" u 2:7 t "check", x/1e6 t"x", b*x**2 t "x^2", f(x) t sprintf("%g + %g * x^2 + %g * x^3", a, b, c)
pause mouse close
EOF
