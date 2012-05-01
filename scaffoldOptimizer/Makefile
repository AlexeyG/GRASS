BNAME = scaffoldOptimizer
OBJ = Configuration.o OverlapperConfiguration.o DPGraph.o DPSolver.o MIQPSolver.o GAIndividual.o GASolver.o FixedMIQPSolver.o ExtendedFixedMIQPSolver.o RelaxedFixedMIQPSolver.o SolverConfiguration.o RandomizedGreedyInitializer.o GAMatrix.o BranchAndBound.o IterativeSolver.o EMSolver.o ScaffoldExtractor.o ScaffoldComparer.o ScaffoldConverter.o GraphViz.o NWAligner.o ContigOverlapper.o optimizer.o
COBJ = Helpers.o DataStore.o DataStoreReader.o Writer.o Timers.o Reader.o ReadCoverageReader.o ReadCoverage.o ReadCoverageRepeatDetecter.o Sequence.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

$(BNAME) : Common core

core : objects
	$(CCC) $(CCCFLAGS) $(CCEFLAGS) $(CCCINC) $(CCEINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB) $(CCELIB)

dirs :
	$(MD) -p $(ODIR)
	$(MD) -p $(BDIR)

objects : $(OBJ)

Common :
	$(MAKE) -C $(CDIR) $(COBJ)

%.o : %.cpp dirs
	$(CCC) -c $(CCCFLAGS) $(CCEFLAGS) $(CCCINC) $(CCEINC) -o $(ODIR)/$@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
