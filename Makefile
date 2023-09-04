P= dk.py
O= out

all:
	python $P > $O.eps
	epstopdf $O.eps
	grep '^%=CSV ' $O.eps | cut -c7- > $O.csv
	grep '^%=OBJ ' $O.eps | cut -c7- > $O.obj
	grep '^%=OFF ' $O.eps | cut -c7- > $O.off
	wc -c -l $O.csv $O.obj $O.off

clean:
	-rm -f $O.*

.PHONY:	all clean

