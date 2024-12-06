-include ../../../petscdir.mk

PETSC_DIR := /home/dai/PETSc/petsc-v3.19.0
PETSC_ARCH := arch-linux-c-debug
# Ensure the PETSC_DIR environment variable is correctly defined.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


# Define compiler flags
CFLAGS := -g -W -Wall -O0 -Wuninitialized
FFLAGS := 
CPPFLAGS := 
FPPFLAGS :=
MESH_DIM ?= 3  # default 3d

# Define source and include directories
CPGRID_PATH := /home/dai/reservoir_simulation/cpgrid/include_src/
MY_SRC_DIR := ${CPGRID_PATH}/src
MY_INCLUDE_DIR := ${CPGRID_PATH}/include

# Define the source files and object files
SRCS := main.C  $(wildcard $(MY_SRC_DIR)/*.C)
OBJS := $(SRCS:.C=.o)

# Define the target binary
TARGET := twophase

# Compilation rule for object files
%.o: %.C
	$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@ -I$(MY_INCLUDE_DIR) -w

# Compilation rule for the main target
$(TARGET): $(OBJS)
	$(CXXLINKER) -o $(TARGET) $^ $(PETSC_LIB)  -I$(MY_INCLUDE_DIR)

# Clean rule to remove generated files
clear:
	rm -rf $(TARGET) $(OBJS) ./Test_output/*.m  *.m *.txt