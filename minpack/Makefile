LIB=minpack
FC=f77
FFLAGS=-O
OBJ = \
	chkder.o        \
	dogleg.o	\
	dpmpar.o	\
	enorm.o	\
	fdjac1.o	\
	fdjac2.o	\
	hybrd.o	\
	hybrd1.o	\
	hybrj.o	\
	hybrj1.o	\
	lmder.o	\
	lmder1.o	\
	lmdif.o	\
	lmdif1.o	\
	lmpar.o	\
	lmstr.o	\
	lmstr1.o	\
	qform.o	\
	qrfac.o	\
	qrsolv.o	\
	r1mpyq.o	\
	r1updt.o	\
	rwupdt.o	\

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

clean:
	rm *.o *.a *.so *.pyc
