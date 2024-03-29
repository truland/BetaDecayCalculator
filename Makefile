#TESTING HOW TO MERGE BRANCHES
#making sure nothing broke

#FOR PROFILING add -g -pg to the compilation flags
#Then run 'mpux ' to profile the functions
#Then run 'gprof mpux > out.data' to get the stats

#FOR THE BRANCH PREDICTIONS add --coverage -ftest-coverage -fbranch-predictions -fprofile-arcs to compilation flags
#then run 'mpux ' to profile the branching
#then run 'gcov sourcefile -m -b > out.data' to get the stats 

# THESE ARE THE DIRECTORIES WHERE THINGS ARE MADE AND BUILT
MKDIR	:= mkdir -p
RMDIR	:= rm -rf

#LAYOUT OF THE BUILD DIRECTORY
INCLUDE	:= ./include
SRC	:= ./src
OBJ	:= ./objs
BIN	:= ./bin
DESTDIR	:= /usr/local/bin/

#REGULAR OPTIONS
CXX		:= g++
OPTIMIZE	:= -O2
WARNINGS	:= -Wall -Wpedantic -Wextra
CPPSTD		:= -std=c++17
CXXFLAGS	+= $(WARNINGS) $(CPPSTD) $(OPTIMIZE) -I$(INCLUDE)

#ROOT CAPABILITY
ROOTCONFIG	:= root-config
CXXFLAGS	+= $(shell $(ROOTCONFIG) --cflags)
LDFLAGS		+= $(shell $(ROOTCONFIG) --ldflags)
LDLIBS		:= $(shell $(ROOTCONFIG) --libs)

#GENERIC BUILD
SRCS	= $(wildcard $(SRC)/*.cpp)
OBJS	= $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(SRCS))
CPA	= $(BIN)/BSCalculator

.PHONY: all clean install

all:	$(CPA)

$(CPA): 	$(OBJS) | $(BIN) Makefile
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS) $(LDFLAGS)

$(OBJ)/%.o:	$(SRC)/%.cpp | $(OBJ) $(INCLUDE)/%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN) $(OBJ):
	$(MKDIR) $@

install: $(CPA)
	sudo install $(CPA) $(DESTDIR)

clean:
	$(RMDIR) $(OBJ) $(BIN)
