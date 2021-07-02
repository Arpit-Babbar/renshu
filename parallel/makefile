CXX       = mpic++
INC_DIR   = ../include
CFLAGS    =  -O3


OBJ = poisson3d.o array2d.o

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

TARGETS = poisson3d

all: $(TARGETS)

#compiling stage(Here)
%.o: $(INC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(CFLAGS) -c $(INC_DIR)/*.cc

#fv2d_var_coeff.o:fv2d_var_coeff.cc array2d.o vtk_anim.o initial_conditions.o
#	$(CXX) $(CFLAGS) -c fv2d_var_coeff.cc

poisson3d: poisson3d.cc array3d.o
	$(CXX) $(CFLAGS) -o $@ $^

clean:
	find . -type f | xargs touch
	rm -f $(TARGETS) *.o

run:
#	$(MAKE)
	mpirun -np 4 ./poisson3d
