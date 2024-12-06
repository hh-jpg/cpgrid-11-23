-include ../../../petscdir.mk

PETSC_DIR := /home/daihuanhe/PETSc/petsc-v3.19.0
PETSC_ARCH := arch-linux-c-debug
# Ensure the PETSC_DIR environment variable is correctly defined.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


# Define compiler flags
CFLAGS := -g -W -Wall -O0 -Wuninitialized 
FFLAGS := 
CPPFLAGS := 
FPPFLAGS :=
# Define the MESH_DIM variable
MESH_DIM ?= 2  # default 3d

# Add the MESH_DIM definition to the compilation flags
CXXFLAGS += -DMESH_DIM=$(MESH_DIM)

# Define source and include directories
CPGRID_PATH := /home/daihuanhe/download/reservoir_simulation/cpgrid/include_src/
MY_SRC_DIR := ${CPGRID_PATH}/src
MY_INCLUDE_DIR := ${CPGRID_PATH}/include 

# Define the source files and object files
SRCS := init.C case_fluid.C res_jacobian.C main.C $(wildcard $(MY_SRC_DIR)/*.C)
OBJS := $(SRCS:.C=.o)

# Define the target binary
TARGET := twophase



# Compilation rule for object files
%.o: %.C
	$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@ -I$(MY_INCLUDE_DIR) -w

# Compilation rule for the main target
$(TARGET): $(OBJS)
	$(CXXLINKER) -o $(TARGET) $^ $(PETSC_LIB)  -I$(MY_INCLUDE_DIR) -I.

# Clean rule to remove generated files
clear:
	rm -rf $(TARGET) $(OBJS)  