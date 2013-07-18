
default : rabbit

rabbit : .FORCE
	$(MAKE) -C src

clean :
	$(MAKE) -C src clean

.FORCE :
