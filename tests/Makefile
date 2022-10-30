ROOT   ?= ../src
CFLAGS ?= -I${HOME}/.local/include -I$(ROOT)
CFLAGS += -Wall -Wextra -Wpedantic -Wredundant-decls -Wpointer-arith -Og -fomit-frame-pointer
LDFLAGS?= -L${HOME}/.local/lib -L./

HEADERS = \
	$(ROOT)/config.h \
	$(ROOT)/types.h \
	$(ROOT)/params.h \
	$(ROOT)/poly.h \
	$(ROOT)/fips202.h \
	$(ROOT)/kem.h \
	$(ROOT)/symmetric.h \
	$(ROOT)/fhe.h \
	$(ROOT)/gpqhe.h
SOURCES = \
	$(ROOT)/types.c \
	$(ROOT)/precomp.c \
	$(ROOT)/reduce.c \
	$(ROOT)/ntt.c \
	$(ROOT)/rns.c \
	$(ROOT)/poly.c \
	$(ROOT)/fips202.c \
	$(ROOT)/canemb.c \
	$(ROOT)/rng.c \
	$(ROOT)/symmetric-shake.c \
	$(ROOT)/sample.c \
	$(ROOT)/he-mem.c \
	$(ROOT)/he-info.c \
	$(ROOT)/he-encode.c \
	$(ROOT)/he-encrypt.c \
	$(ROOT)/he-kem.c \
	$(ROOT)/he-add.c \
	$(ROOT)/he-mult.c \
	$(ROOT)/he-automorphism.c \
	$(ROOT)/he-algo.c

libgpqhe.so: $(HEADERS) $(SOURCES)
#	$(CC) $(CFLAGS) -DRANDOM   $(SOURCES) -fPIC -shared -lgcrypt -lm -o $@
#	$(CC) $(CFLAGS) -DAES256   $(SOURCES) -fPIC -shared -lgcrypt -lm -o $@
	$(CC) $(CFLAGS) -DSUPERCOP $(SOURCES) -fPIC -shared -lgcrypt -lm -o $@

test-precomp: precomp.c libgpqhe.so
	$(CC) $(CFLAGS) precomp.c  $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-ntt: ntt.c libgpqhe.so
	$(CC) $(CFLAGS) ntt.c      $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-crt: crt.c $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -DTEST_CRT -DSUPERCOP $(SOURCES) -fPIC -shared -lgcrypt -lm -o libgpqhe.so
	$(CC) $(CFLAGS) crt.c      $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-polymult: polymult.c libgpqhe.so
	$(CC) $(CFLAGS) polymult.c $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -o $@

test-encode: encode.c libgpqhe.so
	$(CC) $(CFLAGS) encode.c   $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-encdec: encdec.c libgpqhe.so
	$(CC) $(CFLAGS) encdec.c   $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-addmult: addmult.c libgpqhe.so
	$(CC) $(CFLAGS) addmult.c  $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-automorphism: automorphism.c libgpqhe.so
	$(CC) $(CFLAGS) automorphism.c $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-algo: algo.c libgpqhe.so
	$(CC) $(CFLAGS) algo.c         $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

test-misc: misc.c libgpqhe.so
	$(CC) $(CFLAGS) misc.c     $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

libckks.so: $(ROOT)/ntdk.h $(ROOT)/ckks.h $(ROOT)/ntdk.c $(ROOT)/ckks.c
	$(CC) $(CFLAGS) -DSUPERCOP $(ROOT)/ntdk.c $(ROOT)/ckks.c -fPIC -shared -lgcrypt -lm -o $@

test-ckks: ckks_test.c libckks.so
	$(CC) $(CFLAGS) ckks_test.c \
		-DGPQHE_PRECOMP=0 \
		-DGPQHE_DEBUG=0 \
		-DGPQHE_PMU=0 \
		$(LDFLAGS) -lckks -lpmu -lopenblas -llapacke -lgcrypt -lgpg-error -lm -o $@

test-rng: rng.c $(HEADERS) $(SOURCES)
#	$(CC) $(CFLAGS) -DSUPERCOP     $(SOURCES) -fPIC -shared -lgcrypt -lm -o libgpqhe.so
#	$(CC) $(CFLAGS) -DRANDOM       $(SOURCES) -fPIC -shared -lgcrypt -lm -o libgpqhe.so
#	$(CC) $(CFLAGS) rng.c          $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@
	$(CC) $(CFLAGS) -DAES256       $(SOURCES) -fPIC -shared -lgcrypt -lm -o libgpqhe.so
	$(CC) $(CFLAGS) rng.c -DAES256 $(LDFLAGS) -lgpqhe -lpmu -lgcrypt -lm -o $@

clean:
	rm test-*
