ROOT = ..
CFLAGS  += -I$(ROOT)/src -I$(ROOT)/pmu
LDFLAGS ?= -L${ROOT}/lib

test-precomp: precomp.c $(ROOT)/lib/libgpqhe.so
	$(CC) $(CFLAGS) precomp.c  $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-ntt: ntt.c $(ROOT)/lib/libgpqhe.so
	$(CC) $(CFLAGS) ntt.c      $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-crt: crt.c $(ROOT)/src/$(HEADERS) $(ROOT)/src/$(SOURCES)
	cd $(LIBDIR) && mv libgpqhe.so libgpqhe_.so
	cd $(ROOT)/src && $(CC) $(CFLAGS) -DTEST_CRT -DSUPERCOP $(SOURCES) -fPIC -shared -lgcrypt -lm -o $(LIBDIR)/libgpqhe.so
	$(CC) $(CFLAGS) crt.c     $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@
	./$@
	cd $(LIBDIR) && mv libgpqhe_.so libgpqhe.so
	rm $@

test-polymul: polymul.c
	$(CC) $(CFLAGS) polymul.c $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-gpqhe: gpqhe.c
	$(CC) $(CFLAGS) gpqhe.c   $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-algo: algo.c
	$(CC) $(CFLAGS) algo.c    $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-misc: misc.c
	$(CC) $(CFLAGS) misc.c    $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-rng: rng.c
	$(CC) $(CFLAGS) rng.c     $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

include $(ROOT)/src/Makefile

clean:
	rm test-*
