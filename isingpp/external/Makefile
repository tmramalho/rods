CC = g++
SOURCES=$(wildcard *.c)
LIBS=$(SOURCES:.c=.a)
CFLAGS = -Wall

all: ${LIBS}

%.a: %.o
	ar -cr lib$@ $?
	ranlib lib$@
	
%.o: %.c
	${CC} -c ${CFLAGS} -o $@ $?

%.o: %.cpp
	${CC} -c ${CFLAGS} -o $@ $?

clean:
	rm -f *.a *.o *~


