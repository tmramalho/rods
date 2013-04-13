CC = g++
SOURCES=$(wildcard *.cpp)
OBJS=$(SOURCES:.cpp=.o)
CFLAGS = -Wall -g
EXECUTABLE = rods

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): ${OBJS}
	${CC} ${CFLAGS} ${LIBS} -o $@ ${OBJS}

%.o: %.c
	${CC} -c ${CFLAGS} -o $@ $?

%.o: %.cpp
	${CC} -c ${CFLAGS} -o $@ $?

clean:
	rm -f *.o *~ $(EXECUTABLE)

wipe:
	rm -f *.o *~ *.png *.dat *.log $(EXECUTABLE)

