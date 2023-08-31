P= dk.py
O= out

all:	run pdf

run:
	python $P > $O.eps
	grep '^%=CSV ' $O.eps | cut -c7- > $O.csv
	grep '^%=OBJ ' $O.eps | cut -c7- > $O.obj
	grep '^%=OFF ' $O.eps | cut -c7- > $O.off
	wc -c -l $O.csv $O.obj $O.off

pdf:	$O.pdf
	open $O.pdf

$O.pdf:	$O.eps
	epstopdf $O.eps

clean:
	-rm -f $O.*

3:
	cp -p $P 2.py
	cp -p $P 3.py
	2to3- --no-diffs -n -w 3.py
	diff 2.py 3.py | grep -v print

.PHONY:	all run pdf clean 3

