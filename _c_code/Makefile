CC=gcc
CFLAGS= -I. -I/usr/include -L/usr/lib -lgsl -lgslcblas -lm
# run gsl-config --cflags --libs for the correct flags

.PHONY: clean

ella: ella.o
	$(CC) -o ella ella.c $(CFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm ella ella.o
