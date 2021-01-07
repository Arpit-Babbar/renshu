CXX       = g++ #-O3 runs faster.
INC_DIR   =../include
CFLAGS    = -Wall #-O3 Removed optimization to see variables in debugging. Remember to bring it back.


OBJ = fv2d_dirichlet.o array2d.o vtk_anim.o initial_conditions.o

ifeq ($(debug),yes)
	CFLAGS += -DDEBUG
	CFLAGS += -DNDEBUG
	CFLAGS += -Wconversion
	CFLAGS += -Werror
	CFLAGS += -Wextra
	CFLAGS += -pedantic
	CFLAGS += -g
endif

ifeq ($(optimize),yes)
	CXX += -O3
endif

TARGETS = fv2d_dirichlet

all: $(TARGETS)




#compiling stage(Here)
%.o: $(INC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(CFLAGS) -c $(INC_DIR)/*.cc


#fv2d_var_coeff.o:fv2d_var_coeff.cc array2d.o vtk_anim.o initial_conditions.o
#	$(CXX) $(CFLAGS) -c fv2d_var_coeff.cc

fv2d_dirichlet: fv2d_dirichlet.cc array2d.o vtk_anim.o initial_conditions.o 
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	find . -type f | xargs touch
	rm -f $(TARGETS) *.o
	rm -f approximate_solution*.vtk

run:
	./fv2d_dirichlet upwind 0.9 1.0 4 0 
