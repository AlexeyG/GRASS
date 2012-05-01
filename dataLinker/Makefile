BNAME = dataLinker
OBJ = Configuration.o PairedReadConverter.o SequenceConverter.o linker.o
COBJ = Helpers.o DataStore.o Timers.o Reader.o Sequence.o XATag.o DataStoreWriter.o AlignmentReader.o Converter.o Aligner.o AlignerConfiguration.o ReadCoverage.o ReadCoverageWriter.o MummerTilingReader.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

$(BNAME) : Common core

core : objects
	$(CCC) $(CCCFLAGS) $(CCCINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB)

dirs :
	$(MD) -p $(ODIR)
	$(MD) -p $(BDIR)

objects : $(OBJ)

Common :
	$(MAKE) -C $(CDIR) $(COBJ)

%.o : %.cpp dirs
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $(ODIR)/$@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
