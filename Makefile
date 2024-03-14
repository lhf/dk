R= 6
I= -
O= out

all:	dk stats

dk kf:
	python $@.py $R $I > $O.eps
	-pstopdf $O.eps

stats:
	grep '^%=VEF ' $O.eps
	grep '^%=CSV ' $O.eps | cut -c7- > $O.csv
	grep '^%=OBJ ' $O.eps | cut -c7- > $O.obj
	grep '^%=OFF ' $O.eps | cut -c7- > $O.off
	wc -c -l $O.csv $O.obj $O.off

clean:
	-rm -f $O.*

.PHONY:	all dk kf stats clean

