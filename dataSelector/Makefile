BNAME = dataSelector
OBJ = Configuration.o selector.o
COBJ = Helpers.o DataStore.o Timers.o Reader.o Writer.o Sequence.o XATag.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

.PHONY : main Common clean dirs

main : Common $(BDIR)/$(BNAME)

$(BDIR)/$(BNAME) : $(_OBJ) | $(BDIR)
	$(CCC) $(CCCFLAGS) $(CCCINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB)

dirs : | $(ODIR) $(BDIR)

$(ODIR) :
	$(MD) -p $(ODIR)

$(BDIR) :
	$(MD) -p $(BDIR)

Common :
	$(MAKE) -C $(CDIR) $(_COBJ)

$(ODIR)/%.o : %.cpp | $(ODIR)
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
