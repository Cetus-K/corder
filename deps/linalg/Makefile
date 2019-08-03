# Makefile

# compiler
CXX = mpicxx
CFLAGS = -g -O4 -Wall

# file
SRCS = vec.cpp mat.cpp
OBJS = $(SRCS:%.cpp=%.o)

# libraries
LDFLAGS = 
LIBS = 

# includes
INCLUDE = -I$(PREFIX)/deps/includes

# defined macro
.SUFFIXES : .o .cpp

# force compile
all: $(OBJS) liblinalg.a

# compile: -MMD -MP
%.o: %.cpp vec.h mat.h
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ -c $<

liblinalg.a: vec.o mat.o
	ar rsv $@ $?

# clear
clean:
	rm $(OBJS)

.PHONY: all clean

