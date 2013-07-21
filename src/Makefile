
CFLAGS = -Wall -O0
EXE = rabrun rabfile

default : $(EXE)

rabrun : rabrun.o tpl.o

rabfile : rabfile.o rabbit.o tpl.o

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

rabrun.o : rabbit.c
	$(CC) $(CFLAGS) -c -o $@ $<

rabbit.o : rabbit.c
	$(CC) $(CFLAGS) -c -o $@ $< -DRABBIT_LIB

clean :
	$(RM) $(EXE) *.o
