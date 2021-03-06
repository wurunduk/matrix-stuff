CC=g++
CFLAGS=-c -Wall -O3 -ffast-math -march=native
LDFLAGS=-pthread
SOURCES=main.cpp matrix.cpp fpu_init.cpp time.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
