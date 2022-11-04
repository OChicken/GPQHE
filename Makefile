CC		 ?= /usr/bin/gcc
CFLAGS += -Wall -Wextra -Wpedantic -Wredundant-decls -Wpointer-arith -Og -fomit-frame-pointer
LIBDIR ?= $(ROOT)/lib

all:
	cd pmu && $(MAKE) || exit 1
	cd src && $(MAKE) || exit 1
