#!/bin/sh
# Test gribdump
# 

testinput=tests/gribs.wmo

testout=tests/gribs.br
briefout=gribs.br
echo "*** gribdump -b ... \c"
./gribdump -b $testinput > $briefout
if cmp $briefout $testout ; then
    echo "OK"
    rm $briefout
else
    echo "failed!"
    diff $testout $briefout
fi

testout=tests/gribs.dec
decout=gribs.dec
echo "*** gribdump -v ... \c"
./gribdump -v -l /dev/null $testinput > $decout
if cmp $decout $testout ; then
    echo "OK"
    rm $decout
else
    echo "failed!"
    diff $testout $decout
fi

testout=tests/quasi.dec
decout=quasi.dec
echo "*** gribdump -q ... \c"
# to compare floats on different machines, use only 6 significant digits
./gribdump -q "lin,dlat=2.5,dlon=5.0" -v -l /dev/null -p 6 tests/quasi.wmo > $decout
if cmp $decout $testout ; then
    echo "OK"
    rm $decout
else
    echo "failed!"
    diff $testout $decout
fi
