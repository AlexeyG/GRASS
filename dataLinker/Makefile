CCC = g++

CCCFLAGS = -O2 -m64 -std=gnu++0x -Wall
CCCINC = -I../Common/ -I/data/bio/alexeygritsenk/apps/include/
CCCLIB = -L/data/bio/alexeygritsenk/apps/lib/ -lbamtools

MAKE = make
MD = mkdir

CDIR = ../Common
ODIR = obj
BDIR = ../bin

BNAME = dataLinker
OBJ = Configuration.o PairedReadConverter.o SequenceConverter.o linker.o
COBJ = Helpers.o DataStore.o Timers.o Reader.o Sequence.o XATag.o DataStoreWriter.o AlignmentReader.o Converter.o Aligner.o AlignerConfiguration.o ReadCoverage.o ReadCoverageWriter.o MummerTilingReader.o

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

$(BNAME) : Common objects
	$(CCC) $(CCCFLAGS) $(CCCINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB)

dirs :
	$(MD) -p $(ODIR)
	$(MD) -p $(BDIR)

objects : $(OBJ)

Common :
	$(MAKE) -C $(CDIR) -f $(CDIR)/Makefile $(COBJ)

%.o : %.cpp dirs
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $(ODIR)/$@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
