OBJ = Aligner.o AlignmentReader.o DataStore.o DataStoreWriter.o MummerCoordReader.o ReadCoverage.o ReadCoverageRepeatDetecter.o Reader.o Timers.o  XATag.o AlignerConfiguration.o Converter.o DataStoreReader.o Helpers.o MummerTilingReader.o ReadCoverageReader.o ReadCoverageWriter.o Sequence.o Writer.o 

include ../Makefile.config

Common : objects

dirs :
	$(MD) -p $(ODIR)

objects : $(OBJ)

%.o : %.cpp dirs
	$(CCC) -c $(CCCFLAGS) $(CCCINC) -o $(ODIR)/$@ $<

clean :
	$(RM) -rf $(ODIR)
