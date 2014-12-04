
CC = gcc
CFLAGS = -O3
LIBS = -lm

HEADERS = func.h
OBJECTS = H2_main.o func.o
PROGRAM = H2

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

