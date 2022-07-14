include ../make.inc

all:  testlsame testslamch testdlamch testsecond testdsecnd 

testlsame:  lsame.o lsametst.o
	$(LOADER) -o testlsame lsame.o lsametst.o

testslamch: slamch.o lsame.o slamchtst.o
	$(LOADER) -o testslamch slamch.o lsame.o slamchtst.o

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(LOADER) -o testdlamch dlamch.o lsame.o dlamchtst.o

testsecond: second.o secondtst.o
	$(LOADER) -o testsecond second.o secondtst.o

testdsecnd: dsecnd.o dsecndtst.o
	$(LOADER) -o testdsecnd dsecnd.o dsecndtst.o

slamch.o: slamch.f ; $(FORTRAN) $(NOOPT) -c $<
dlamch.o: dlamch.f ; $(FORTRAN) $(NOOPT) -c $<

.f.o: ; $(FORTRAN) $(OPTS) -c $<
