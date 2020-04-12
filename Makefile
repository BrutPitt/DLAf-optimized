CC = g++
TARGET = dlaf
COMPILE_FLAGS = -std=c++14 -flto -O3 -Wall -Wextra -pedantic -Wno-unused-parameter -march=native

all: $(TARGET)

$(TARGET): $(TARGET).o
	$(CC) -o $(TARGET) $(TARGET).o 

$(TARGET).o: $(TARGET).cpp
	$(CC) -c $(COMPILE_FLAGS) $(TARGET).cpp

clean:
	$(RM) $(TARGET)
