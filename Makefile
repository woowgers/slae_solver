CC = gcc
CFLAGS = -g3 -lm

.PHONY: main test clean

all: main limits rational test


main: src/main.c
	$(CC) -o main src/slae.c src/rational.c src/main.c $(CFLAGS)

limits: src/limits.c
	$(CC) -o limits src/limits.c $(CFLAGS)

rational: src/test_rational.c
	$(CC) -o rational src/test_rational.c src/rational.c $(CFLAGS)

build_tests: src/test.c
	$(CC) -o test src/slae.c src/test.c -lcriterion

test: build_tests
	@./test

clean:
	@rm -f rational limits main

