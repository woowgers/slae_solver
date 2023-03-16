CC = gcc
CFLAGS = -g3 -lm

.PHONY: main test clean

all: main limits rational test_remove test stepik


main: src/main.c
	$(CC) -o bin/main src/slae.c src/rational.c src/main.c $(CFLAGS)

limits: src/limits.c
	$(CC) -o bin/limits src/limits.c $(CFLAGS)

rational: src/test_rational.c
	$(CC) -o bin/rational src/test_rational.c src/rational.c $(CFLAGS)

test_remove: src/test_remove.c
	$(CC) -o bin/test_remove src/test_remove.c src/slae.c $(CLFAGS)


stepik: src/slae.h src/slae.h build_stepik/stepik.def.c
	cat src/slae.h src/slae.c build_stepik/stepik.def.c > build_stepik/stepik.c
	sed -i 's/#include ".*//g' build_stepik/stepik.c
	$(CC) -o build_stepik/stepik build_stepik/stepik.c $(CFLAGS)


build_tests: src/test.c
	$(CC) -o bin/test src/slae.c src/test.c -lcriterion $(CFLAGS)

test: build_tests
	bin/test


clean:
	rm -f bin/*
	rm -f build_stepik/stepik.c
