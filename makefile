CXX      = g++
CXXFLAGS = -Wall -march=native -O3 -flto -fuse-linker-plugin -fopenmp
LDFLAGS  = -lgsl -lgslcblas

OBJECTS = main.o settings.o lattice.o scc_calc.o plot.o
DEFINES = -D_EIGEN_DONT_PARALLELIZE

mfhub : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o mfhub

main.o : main.cpp typedefs.hpp settings.hpp scc_inout.hpp scc_calc.hpp plot.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o

settings.o : settings.hpp settings.cpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c settings.cpp -o settings.o
	
lattice.o : lattice.hpp lattice.cpp typedefs.hpp settings.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c lattice.cpp -o lattice.o
	
scc_calc.o : scc_calc.hpp scc_calc.cpp typedefs.hpp settings.hpp scc_inout.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c scc_calc.cpp -o scc_calc.o
	
plot.o : plot.hpp plot.cpp typedefs.hpp settings.hpp scc_inout.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c plot.cpp -o plot.o

clean:
	rm mfhub $(OBJECTS)
