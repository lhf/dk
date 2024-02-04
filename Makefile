P= dk.py
O= out
R= 6
I= -

all:
	python $P $R $I > $O.eps
	grep '^%=CSV ' $O.eps | cut -c7- > $O.csv
	grep '^%=OBJ ' $O.eps | cut -c7- > $O.obj
	grep '^%=OFF ' $O.eps | cut -c7- > $O.off
	wc -c -l $O.csv $O.obj $O.off
	-pstopdf $O.eps

clean:
	-rm -f $O.*

.PHONY:	all clean

