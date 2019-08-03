# Makefile

# compiler
CXX = mpicxx
CFLAGS = -g -O4 -Wall

# file
SRCS = function.cpp
OBJS = $(SRCS:%.cpp=%.o)

# libraries
LDFLAGS = 
LIBS = 

# includes
INCLUDE = 

# defined macro
.SUFFIXES : .o .cpp

# force compile
all: $(OBJS) libfunction.a

# compile: -MMD -MP
%.o: %.cpp header.h
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ -c $<

libfunction.a: function.o
	ar rsv $@ $?

# clear
clean:
	rm $(OBJS)

.PHONY: all clean

