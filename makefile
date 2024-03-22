CXX := g++
CXXFLAGS := -O3 -std=c++11 -Wall -Wconversion -Wextra -Wpedantic

TARGET := main
OBJS := main.o CGSolver.o matvecops.o COO2CSR.o heat.o sparse.o
INCS := CGSolver.hpp matvecops.hpp COO2CSR.hpp heat.hpp sparse.hpp

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS)

%.o: %.cpp $(INCS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) *~
