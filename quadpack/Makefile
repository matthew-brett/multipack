LIB=quadpack
FFLAGS=-O
OBJ=\
dqag.o\
dqage.o\
dqagi.o\
dqagie.o\
dqagp.o\
dqagpe.o\
dqags.o\
dqagse.o\
dqawc.o\
dqawce.o\
dqawf.o\
dqawfe.o\
dqawo.o\
dqawoe.o\
dqaws.o\
dqawse.o\
dqc25c.o\
dqc25f.o\
dqc25s.o\
dqcheb.o\
dqelg.o\
dqk15.o\
dqk15i.o\
dqk15w.o\
dqk21.o\
dqk31.o\
dqk41.o\
dqk51.o\
dqk61.o\
dqmomo.o\
dqng.o\
dqpsrt.o\
dqwgtc.o\
dqwgtf.o\
dqwgts.o\

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a 

install:	lib$(LIB).a
	mv lib$(LIB).a /usr/local/lib
	rm *.o
