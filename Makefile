SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#---------------------------------------------------------------------------------------------------
#
# Set CPLEXDIR and CONCERTDIR to the directories where CPLEX and CONCERT are installed.
#
#---------------------------------------------------------------------------------------------------

CPLEXDIR      = <local_cplex_path>
CONCERTDIR    = <local_concert_path>

#---------------------------------------------------------------------------------------------------
# Compiler selection
#---------------------------------------------------------------------------------------------------

CXX = g++

#---------------------------------------------------------------------------------------------------
# Directories
#---------------------------------------------------------------------------------------------------

OBJDIR = build
SRCDIR = src

#---------------------------------------------------------------------------------------------------
# Executables
#---------------------------------------------------------------------------------------------------

EXE = SepComp S_MIP CombClustOut

#---------------------------------------------------------------------------------------------------
# Object files
#---------------------------------------------------------------------------------------------------
SEPOBJ = SeperateGraphIntoComponents.o Graph.o Vertex.o Edge.o
SOBJ = S_MIP.o Graph.o Vertex.o Edge.o Timer.o
CCOOBJ = CombClustOut.o Graph.o Vertex.o Edge.o

#---------------------------------------------------------------------------------------------------
# Compiler options
#---------------------------------------------------------------------------------------------------

CXXFLAGS = -O3 -fPIC -fexceptions -DIL_STD -std=c++11#-fno-strict-aliasing

#---------------------------------------------------------------------------------------------------
# Link options and libraries
#---------------------------------------------------------------------------------------------------

CPLEXLIBDIR    = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR  = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CXXLNDIRS      = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CXXLNFLAGS     = -lconcert -lilocplex -lcplex -lm -lpthread -ldl

#---------------------------------------------------------------------------------------------------
# Includes
#---------------------------------------------------------------------------------------------------

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

INCLUDES = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

#---------------------------------------------------------------------------------------------------
all: CXXFLAGS += -DNDEBUG
all: $(EXE)

debug: CXXFLAGS += -g
debug: $(EXE)

SepComp: $(addprefix $(OBJDIR)/, SeperateGraphIntoComponents.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(SEPOBJ))

S_MIP: $(OBJDIR)/S_MIP.o
	$(CXX) $(CXXLNDIRS) -o $@ $(addprefix $(OBJDIR)/, $(SOBJ)) $(CXXLNFLAGS)

CombClustOut: $(OBJDIR)/CombClustOut.o
	$(CXX) $(CXXLNDIRS) -o $@ $(addprefix $(OBJDIR)/, $(CCOOBJ)) $(CXXLNFLAGS)

$(OBJDIR)/SeperateGraphIntoComponents.o: $(addprefix $(SRCDIR)/, SeperateGraphIntoComponents.cpp) \
																																									$(addprefix $(OBJDIR)/, Graph.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/S_MIP.o: $(addprefix $(SRCDIR)/, S_MIP.cpp) \
																			$(addprefix $(OBJDIR)/, Graph.o Vertex.o Edge.o Timer.o)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)/CombClustOut.o: $(addprefix $(SRCDIR)/, CombClustOut.cpp)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Timer.o: $(addprefix $(SRCDIR)/, Timer.cpp Timer.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Graph.o: $(addprefix $(SRCDIR)/, Graph.cpp Graph.h) \
																			$(addprefix $(OBJDIR)/, Vertex.o Edge.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Vertex.o:	$(addprefix $(SRCDIR)/, Vertex.cpp Vertex.h ) \
																				$(addprefix $(OBJDIR)/, Edge.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Edge.o:	$(addprefix $(SRCDIR)/, Edge.cpp Edge.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#---------------------------------------------------------------------------------------------------
.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o
#---------------------------------------------------------------------------------------------------
