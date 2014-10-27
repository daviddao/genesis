# --------------------------------
#   Program
# --------------------------------

PROGRAM = genesis

# --------------------------------
#   Compiler Options
# --------------------------------

# Debug & Profiling (comment out if not needed)
DBG     = -g -pg -DDEBUG

# Warning flags
WARN    = -Wall -Wextra -pedantic-errors

# Compiler flags
STDCC   = g++
MPICC   = mpic++
CCFLAGS = -std=c++11 -O1 $(WARN) $(DBG)
LDFLAGS = -lm
#-lpll-sse3 -lm

# --------------------------------
#   File lists
# --------------------------------

SRCFILES := $(shell find ./src -type f -name "*.cc")
HDRFILES := $(shell find ./src -type f -name "*.hh")
OBJFILES := $(patsubst %.cc,%.o,$(SRCFILES))
DEPFILES := $(patsubst %.cc,%.d,$(SRCFILES))

ALLFILES := $(SRCFILES) $(HDRFILES)

# --------------------------------
#   Make rules
# --------------------------------

.PHONY: all clean dist check test todo

# Include dependecies
# (they are generated when compiling the sources and contain makefile-formatted
# information on their dependecies, so that nothing should be missing in the end)
-include $(DEPFILES)

# Build the standard version of the program
all: CC = ${STDCC}
all: $(PROGRAM)
	@echo "\n========== Done std  =========="

# Build an MPI version of the program
mpi: CC = ${MPICC}
mpi: CCFLAGS += -DUSE_MPI
mpi: $(PROGRAM)
	@echo "\n========== Done mpi  =========="

# Link all objects to get the program
$(PROGRAM): $(OBJFILES)
	@echo "\n==========  Linking  =========="
	@echo "Objects: $(OBJFILES)\n"
	$(CC) $(OBJFILES) -o $@ $(LDFLAGS)

# Compile all sources
%.o: %.cc Makefile
	@echo "\n========== Compiling =========="
	@echo "File: $< > $@"
	$(CC) $(CCFLAGS) -MMD -MP -c $< -o $@

# Remove all generated files
clean:
	@echo "\n========== Cleaning  =========="
	-@$(RM) $(PROGRAM) $(OBJFILES) $(DEPFILES)

# Extract todos
todo:
	-@$(RM) TODO
	-@for file in $(ALLFILES:Makefile=); do fgrep -HnT -e TODO -e FIXME $$file \
	| sed "s/[[:space:]]*[\/\*]*[[:space:]]*TODO[[:space:]]*/ /g" >> TODO; done; true
