
CFLAGS = -Wall -O3
EXE = rabrun rabfile rabview hydro hydro-2d

GL_LIB = -framework OpenGL -framework GLUT
GL_INC = -I/usr/X11R6/include

default : $(EXE)

rabrun : rabrun.o tpl.o
	$(CC) $(CFLAGS) -o $@ $^

rabfile : rabfile.o rabbit.o tpl.o
	$(CC) $(CFLAGS) -o $@ $^

rabview : rabview.o rabbit.o tpl.o
	$(CC) $(CFLAGS) -o $@ $^ $(GL_LIB)

hydro : hydro.o rabbit.o tpl.o euler1d.o
	$(CC) $(CFLAGS) -o $@ $^

hydro-2d : hydro-2d.o rabbit.o tpl.o euler1d.o
	$(CC) $(CFLAGS) -o $@ $^

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

rabview.o : rabview.c
	$(CC) $(CFLAGS) -c -o $@ $< $(GL_INC)

rabrun.o : rabbit.c
	$(CC) $(CFLAGS) -c -o $@ $<

rabbit.o : rabbit.c
	$(CC) $(CFLAGS) -c -o $@ $< -DRABBIT_LIB

clean :
	$(RM) $(EXE) *.o
