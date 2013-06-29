
EXE = ztree-test
CFLAGS = -Wall -g -O0
CC = cc

default : $(EXE)

ztree-test : ztree-test.o ztree.o zmesh.o

clean :
	$(RM) *.o $(EXE)
