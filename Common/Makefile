OBJ = Aligner.o AlignmentReader.o DataStore.o DataStoreWriter.o MummerCoordReader.o ReadCoverage.o ReadCoverageRepeatDetecter.o Reader.o Timers.o  XATag.o AlignerConfiguration.o Converter.o DataStoreReader.o Helpers.o MummerTilingReader.o ReadCoverageReader.o ReadCoverageWriter.o Sequence.o Writer.o 

include ../Makefile.config

_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))

.PHONY : objects clean dirs

objects : $(_OBJ)

dirs : | $(ODIR)

$(ODIR) :
	$(MD) -p $(ODIR)

$(ODIR)/%.o : %.cpp %.h | $(ODIR) 
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $@ $<

clean :
	$(RM) -rf $(ODIR)
