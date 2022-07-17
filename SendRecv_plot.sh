#!/bin/bash

if [ $# != 1 ] ; then
	echo "Usage: $0 file_from_fatring"
	exit 1
fi

gnuplot << EOF
se log xy
se xla "#chunks"
se yla "time (s)"
p [][1e-5:100] "$1" u 2:4 t "setup", "" u 2:5 t "SendRecv", "" u 2:6 t "Waitall", "" u 2:7 t "check", x/1e6 t"x", x**2/1e9 t "x^2"
pause mouse close
EOF
