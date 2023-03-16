CC = gcc
CFLAGS = -g3 -lm

.PHONY: main test clean

all: main limits rational test


main: src/main.c
	$(CC) -o bin/main src/slae.c src/rational.c src/main.c $(CFLAGS)

limits: src/limits.c
	$(CC) -o bin/limits src/limits.c $(CFLAGS)

rational: src/test_rational.c
	$(CC) -o bin/rational src/test_rational.c src/rational.c $(CFLAGS)

build_tests: src/test.c
	$(CC) -o bin/test src/slae.c src/test.c -lcriterion

test: build_tests
	@bin/test

clean:
	@rm -f rational limits main

