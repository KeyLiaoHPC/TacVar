# Build FilT
CC ?= gcc
CFLAGS ?= -O2 -std=c99 -Wall
LDFLAGS ?= -lm

.PHONY: all clean

all: filt.x

filt.x: filt.c filt.h
	$(CC) $(CFLAGS) -o $@ filt.c $(LDFLAGS)

clean:
	rm -f filt.x
