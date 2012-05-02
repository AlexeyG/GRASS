BNAME = dataLinker
OBJ = Configuration.o PairedReadConverter.o SequenceConverter.o linker.o
COBJ = Helpers.o DataStore.o Timers.o Reader.o Sequence.o XATag.o DataStoreWriter.o AlignmentReader.o Converter.o Aligner.o AlignerConfiguration.o ReadCoverage.o ReadCoverageWriter.o MummerTilingReader.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

.PHONY : main Common objects clean dirs

main : Common $(BDIR)/$(BNAME)

$(BDIR)/$(BNAME) : objects $(BDIR)
	$(CCC) $(CCCFLAGS) $(CCCINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB)

dirs : $(ODIR) $(BDIR)

$(ODIR) :
	$(MD) -p $(ODIR)

$(BDIR) :
	$(MD) -p $(BDIR)

objects : $(_OBJ)

Common :
	$(MAKE) -C $(CDIR) $(_COBJ)

$(ODIR)/%.o : %.cpp $(ODIR)
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
