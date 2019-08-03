# Makefile

# compiler
CXX = mpicxx
CFLAGS = -g -O4 -Wall
INCLUDES_DIR = $(PREFIX)/deps/includes
LINALG_DIR = $(PREFIX)/deps/linalg
VORO_DIR = $(PREFIX)/deps/voro++
QHULL_DIR = $(PREFIX)/deps/qhull

# file
SRCS = corder.cpp param.cpp gofr.cpp sofq.cpp dfc.cpp covoro.cpp lboo.cpp bacf.cpp pofqw.cpp bofree.cpp povoro.cpp coqhull.cpp csform.cpp
OBJS = $(SRCS:%.cpp=%.o)
TARGET = corder

# libraries
LDFLAGS = -L$(INCLUDES_DIR) -L$(LINALG_DIR) -L$(VORO_DIR)/lib -L$(QHULL_DIR)/lib
LIBS = -lfunction -llinalg -lvoro++ -lqhullcpp -lqhullstatic_r -lirc

# includes
INCLUDE = -I$(INCLUDES_DIR) -I$(LINALG_DIR) -I$(VORO_DIR)/src -I$(QHULL_DIR)/src

# defined macro
.SUFFIXES : .o .cpp

# force compile
all: $(TARGET)

# link
$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o $@ $^ $(LIBS)

# compile: -MMD -MP
%.o: %.cpp corder.h
	$(CXX) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o $@ -c $< $(LIBS)

# clear
clean:
	rm $(TARGET) $(OBJS)

.PHONY: all clean

