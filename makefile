TARGET = exe
CC = g++
LIBS = -lm -lboost_filesystem -lboost_iostreams -lboost_system
HEAD = ./include
SRCS = ./source
INCDIR = -I$(HEAD)
CFLAGS = -Wall -ggdb3 -std=c++17 -funroll-loops $(INCDIR)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp))

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	-rm -r $(SRCS)/*.o
	-rm -r $(TARGET)

