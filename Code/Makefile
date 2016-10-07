OBJS =  main.o krig_utility.o model.o lhs.o chol.o chol_sol.o cal_mean.o predictor.o likelihood.o mads_poll.o msecheck.o push_to_grid.o cmaes_driver.o 

# Compilers 

CC = gcc
CXX = g++

# Flags

CXXFLAGS = -std=c++11
LDFLAGS = -O3  

# Commands

COMMAND = krig

%.o: %.cpp 
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

$(COMMAND): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(COMMAND) $(OBJS)

clean: 
	rm -f $(OBJS) $(COMMAND)
