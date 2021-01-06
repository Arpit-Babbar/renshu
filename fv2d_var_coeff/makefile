CXX       = g++ #-O3 runs faster.
INC_DIR   =../include
CFLAGS    = -Wall #-O3 Removed optimization to see variables in debugging. Remember to bring it back.


OBJ = fv2d_var_coeff.o array2d.o vtk_anim.o initial_conditions.o

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

TARGETS = fv2d_var_coeff

all: $(TARGETS)

#array2d.o: $(INC_DIR)/array2d.h $(INC_DIR)/array2d.cc 
#initial_conditions.o: $(INC_DIR)/initial_conditions.h $(INC_DIR)/initial_conditions.cc 
#vtk_anim.o: $(INC_DIR)/vtk_anim.h $(INC_DIR)/vtk_anim.cc 



#compiling stage(Here)
%.o: $(INC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(CFLAGS) -c $(INC_DIR)/*.cc


#fv2d_var_coeff.o:fv2d_var_coeff.cc array2d.o vtk_anim.o initial_conditions.o
#	$(CXX) $(CFLAGS) -c fv2d_var_coeff.cc

fv2d_var_coeff: fv2d_var_coeff.cc array2d.o vtk_anim.o initial_conditions.o 
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	find . -type f | xargs touch
	rm -f $(TARGETS) *.o
	rm -f approximate_solution*.vtk

run:
	./fv2d_var_coeff lw 0.9 2pi 4 0 
