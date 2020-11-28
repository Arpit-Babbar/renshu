CXX = g++ -Wall

TARGETS = limiter

all: $(TARGETS)

vector_upgrade.o: vector_upgrade.cc vector_upgrade.h
initial_conditions.o: initial_conditions.cc initial_conditions.h
finite_volume_solver.o: finite_volume_solver.cc finite_volume_solver.h

limiter: limiter.o vector_upgrade.o initial_conditions.o finite_volume_solver.o #making a file without extension?
	$(CXX) -o $@ $^

#linsol_test: sparse_matrix.o Vector.o cg_solver.o jacobi_solver.o \
#             sor_solver.o ssor_solver.o linsol_test.o
#	$(CXX) -o $@ $^

#linsol_test_2d: sparse_matrix.o Vector.o cg_solver.o jacobi_solver.o \
                sor_solver.o ssor_solver.o linsol_test_2d.o
#	$(CXX) -o $@ $^

clean:
	rm -f $(TARGETS) *.o

run:
	/limiter
