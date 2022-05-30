#=============================================================================
TOP= ./
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/
OUT= $(TOP)output/

SYSLIB= -lm

VPATH= $(SRC)
#=============================================================================
CC=g++#icpc#
CFLAGS= -g -O2
#==========================================================================
ifeq ($(CC),g++)
	CFLAGS+= -std=c++14 -Wall -Wextra #-fopenmp
		#-fcheck=all
endif

ifeq ($(CC),icpc)
	CFLAGS+= -std=c++14 -Wall
		#-check-bounds
endif
#=============================================================================
OBJ= $(addprefix $(OBJDIR), \
	main.o \
	grid_data.o \
	field.o \
	solve_metric_fields.o \
	initial_data.o \
	diagnostics.o \
	evolve_scalar_field.o \
	outputfiles.o \
	fd_stencils.o	\
	compute_potentials.o \
	)
DEPS=		grid_data.hpp	field.hpp solve_metric_fields.hpp	initial_data.hpp	diagnostics.hpp evolve_scalar_field.hpp	outputfiles.hpp	fd_stencils.hpp compute_potentials.hpp
#=============================================================================
all: default.run
test: $(TEST)
#=============================================================================
%.run: $(OBJ)
	$(CC) $^ -o $(BIN)$@ $(SYSLIB) $(CFLAGS)
#=============================================================================
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<
#=============================================================================
# .PHONY: clean clean_out
# clean:
# 	$(RM) $(OBJDIR)*.o
# 	$(RM) $(BIN)*.run
# clean_out:
# 	$(RM) -r $(OUT)*
