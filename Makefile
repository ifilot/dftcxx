# set compiler and compile options
EXEC = dftcxx
CXX = g++                             # use the GNU C++ compiler
OPTS = -std=c++0x -O3 -Wall           # use some optimization, report all warnings and enable debugging
CFLAGS = $(OPTS)                      # add compile flags
LDFLAGS = # specify link flags here

# set a list of directories
INCDIR = ./include
OBJDIR = ./obj
BINDIR = ./bin
SRCDIR = ./src

# set the include folder where the .h files reside
CFLAGS += -I$(INCDIR) -I$(SRCDIR)

# specify OS-related stuff
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CFLAGS  += -I/opt/local/include/eigen3 -I/opt/local/include
	LDFLAGS +=
else
    CFLAGS += `pkg-config --cflags eigen3`
	LDFLAGS +=
endif

# add here the source files for the compilation
SOURCES = dftcxx.cpp cgf.cpp integrals.cpp gamma.cpp quadrature.cpp \
molecule.cpp moleculargrid.cpp dft.cpp functionals.cpp

# create the obj variable by substituting the extension of the sources
# and adding a path
_OBJ = $(SOURCES:.cpp=.o)
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

all: $(BINDIR)/$(EXEC)

$(BINDIR)/$(EXEC): $(OBJ)
	$(CXX) -o $(BINDIR)/$(EXEC) $(OBJ) $(LIBDIR) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c -o $@ $< $(CFLAGS)

clean:
	rm -vf $(BINDIR)/$(EXEC) $(OBJ)
