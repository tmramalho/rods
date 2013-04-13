CC = g++
SOURCES=len.cpp sim.cpp lab.cpp main.cpp
OBJS=$(SOURCES:.cpp=.o)
INCDIR = $(CURDIR)/external
LIBDIR = $(CURDIR)/external
LIBS = -lmt19937ar
CFLAGS = -Wall -g
EXECUTABLE = rods

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): ${OBJS}
	${CC} ${CFLAGS} -I${INCDIR} -L${LIBDIR} -o $@ $? ${LIBS}

%.o: %.c
	${CC} -c ${CFLAGS} -I${INCDIR} -o $@ $?

%.o: %.cpp
	${CC} -c ${CFLAGS} -I${INCDIR} -o $@ $?

clean:
	rm -f *.o *~ $(EXECUTABLE)

wipe:
	rm -f *.o *~ *.png *.dat *.log $(EXECUTABLE)

