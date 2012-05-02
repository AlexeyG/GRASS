BNAME = readCleaner
OBJ = Configuration.o PairedReadProcessor.o cleaner.o
COBJ = Helpers.o DataStore.o Timers.o Reader.o Writer.o Sequence.o XATag.o DataStoreWriter.o AlignmentReader.o Converter.o Aligner.o AlignerConfiguration.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

.PHONY : main Common objects clean dirs

main : Common $(BDIR)/$(BNAME)

$(BDIR)/$(BNAME) : objects
	$(CCC) $(CCCFLAGS) $(CCCINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB)

dirs : $(ODIR) $(BDIR)

$(ODIR) :
	$(MD) -p $(ODIR)

$(BDIR) :
	$(MD) -p $(BDIR)

objects : $(_OBJ)

Common :
	$(MAKE) -C $(CDIR) $(_COBJ)

$(ODIR)/%.o : %.cpp dirs
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
