# This is a more elaborate GNU Makefile using:
# - implicit rules,
# - automatic variables.

N ?= 1

# This is a good practice to store values subject to change in a variable.
# Variable can be assigned simply as below and their value returned with $(...)
CC=g++
CFLAGS="-O2"
TARGET=toto
LIB=libq1.a
MAINSRC=main.cpp

# A list of source file can be conveniently created with a wildcard which will
# match all the file with a .c extension. 
#SRC=$(wildcard *.cpp)

SRC=zeta0.cpp mach0.cpp

# As objects do not exist yet, the list can be created from the source list by
# substituting the .c extension with .o
OBJ=$(patsubst %.cpp, %.o, $(SRC))

# The object file for the main source is obtained by changing the extension
MAINOBJ=$(patsubst %.cpp, %.o, $(MAINSRC))
# The list of objects in the library is constructed by omitting main object
LIBOBJ=$(patsubst $(MAINOBJ),,$(OBJ))

# This implicit rule will match any dependency ending with .o
%.o:%.c
	$(CC) $(CFLAGS) -c $<

# The target creates the executable: it depends on a library and a main object
$(TARGET): $(LIB) $(MAINOBJ)
	$(CC) -o $(TARGET) $(MAINOBJ) $(LIB)

all: $(TARGET)

# The library target depends on a list of objects, calling the previous rule
$(LIB): $(LIBOBJ)
	@ar r $@ $(LIBOBJ)

# Another target to cleanup files
clean:
	rm -f $(TARGET) $(LIB) $(OBJ)

# List objects in the library
list_libobjs: $(LIB)
	@echo "All the objects contained in the static library are printed below:\n"
	nm $(LIB)

# Run test
test: $(TARGET)
	./$(TARGET) 0 $(N)

# Run unit test
utest: $(TARGET) 
	./$(TARGET) 1 $(N) 

# Run verification test
vtest: $(TARGET) 
	./$(TARGET) 2 $(N)


