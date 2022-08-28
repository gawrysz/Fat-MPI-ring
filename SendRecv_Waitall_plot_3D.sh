#!/bin/bash

if [ $# != 1 ] ; then
	echo "Usage: $0 file_from_fatring"
	exit 1
fi

gnuplot << EOF
se log
se xla "#chunks"
se yla "#doubles"
se zla "time (s)"

a=1
b=1e-10
c=1e-10
f(x, y) = a +b*y +c*x**2.5
fit              f(x, y) "$1" u 2:(\$2+\$3):6 via a, b, c
fit [:100]       f(x, y) "$1" u 2:(\$2*\$3):6 via a, b
fit [:100][:100] f(x, y) "$1" u 2:(\$2*\$3):6 via a

se tit sprintf("Waitall execution time\nFit: %g + %g * #doubles + %g * #chunks^{2.5}", a, b, c)

sp "$1" u 2:(\$2*\$3):6, f(x,y)
pause mouse close
EOF
