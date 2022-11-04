ROOT ?= ..

libpmu.so: pmu.c
	$(CC) $(CFLAGS) -fPIC -shared $^ -o $(LIBDIR)/$@

include $(ROOT)/Makefile
