#***********************************************#
# Define macro "PRINTING" as follows:		#
# NONE    - printing is disabled		#
# COMPACT - for basic information printing	#
# FULL    - thorough printing			#
#***********************************************#
CC	= g++

DATE	= `date +%d%m%Y-%H%M`

#MACROFLAGS += -DNONE=0 -DCOMPACT=1 -DFULL=2
#MACROFLAGS += -DPRINTING=FULL

ifeq ($(CC),icc)
	CFLAGS	= -O3
else ifeq ($(CC),g++)
	CFLAGS	= -O3 -Wno-unused-result
endif

#CFLAGS	+= $(MACROFLAGS)
#CFLAGS	+= -fopenmp
#CFLAGS	+= -g
#CFLAGS	+= -std=c++0x
#CFLAGS	+= -static
#CFLAGS	+= -pg
#CFLAGS	+= -Wall

ROOT		= $(CURDIR)
SRCDIR		= $(ROOT)/src
HPPDIR		= $(ROOT)/hpp
OBJDIR		= $(ROOT)/obj

CHECK_CC := $(shell which $(CC))
ifeq ($(CHECK_CC),)
	$(error No $(CC) found!)
endif

.PHONY: all clean dirs

all: dirs SpLap

dirs:
	mkdir -p $(OBJDIR)

SpLap: $(OBJDIR)/graph.o $(OBJDIR)/hist.o $(OBJDIR)/load_input.o $(OBJDIR)/main.o $(OBJDIR)/rand_gen.o $(OBJDIR)/Sparsify.o $(OBJDIR)/aux.o $(OBJDIR)/ApproxReff.o
	$(CC) $(OBJDIR)/graph.o $(OBJDIR)/hist.o $(OBJDIR)/load_input.o $(OBJDIR)/main.o $(OBJDIR)/rand_gen.o $(OBJDIR)/Sparsify.o $(OBJDIR)/aux.o $(OBJDIR)/ApproxReff.o -o SpLap

$(OBJDIR)/load_input.o: $(SRCDIR)/load_input.cpp
	$(CC) -c $(SRCDIR)/load_input.cpp $(CFLAGS) -o $(OBJDIR)/load_input.o 
	
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(CC) -c $(SRCDIR)/main.cpp $(CFLAGS) -o $(OBJDIR)/main.o 

$(OBJDIR)/hist.o: $(SRCDIR)/hist.cpp
	$(CC) -c $(SRCDIR)/hist.cpp $(CFLAGS) -o $(OBJDIR)/hist.o

$(OBJDIR)/graph.o: $(SRCDIR)/graph.cpp
	$(CC) -c $(SRCDIR)/graph.cpp $(CFLAGS) -o $(OBJDIR)/graph.o

$(OBJDIR)/rand_gen.o: $(SRCDIR)/rand_gen.cpp
	$(CC) -c $(SRCDIR)/rand_gen.cpp $(CFLAGS) -o $(OBJDIR)/rand_gen.o

$(OBJDIR)/Sparsify.o: $(SRCDIR)/Sparsify.cpp
	$(CC) -c $(SRCDIR)/Sparsify.cpp $(CFLAGS) -o $(OBJDIR)/Sparsify.o

$(OBJDIR)/aux.o: $(SRCDIR)/aux.cpp
	$(CC) -c $(SRCDIR)/aux.cpp $(CFLAGS) -o $(OBJDIR)/aux.o

$(OBJDIR)/ApproxReff.o: $(SRCDIR)/ApproxReff.cpp
	$(CC) -c $(SRCDIR)/ApproxReff.cpp $(CFLAGS) -o $(OBJDIR)/ApproxReff.o

clean:
	rm -rf *~ $(SRCDIR)/*~ $(OBJDIR) SpLap
