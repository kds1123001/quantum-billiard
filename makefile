FC = gfortran
FFLAGS = -O3 -J build
LIBS = -llapack -lblas
SRC = src/kinds.f90 src/geometry.f90 src/hamiltonian.f90 src/evolve.f90 src/render.f90 src/main.f90
BIN = qbilliard_term

all: $(BIN)

$(BIN): $(SRC)
	mkdir -p build
	$(FC) $(FFLAGS) -o $(BIN) $(SRC) $(LIBS)

run: $(BIN)
	./$(BIN)

clean:
	rm -rf build $(BIN)
