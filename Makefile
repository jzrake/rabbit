
EXE = ztree-test

default : $(EXE)

ztree-test : ztree-test.o ztree.o

clean :
	$(RM) *.o $(EXE)
