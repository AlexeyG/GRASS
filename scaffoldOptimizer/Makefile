BNAME = scaffoldOptimizer
OBJ = Configuration.o OverlapperConfiguration.o DPGraph.o DPSolver.o MIQPSolver.o GAIndividual.o GASolver.o FixedMIQPSolver.o ExtendedFixedMIQPSolver.o RelaxedFixedMIQPSolver.o SolverConfiguration.o RandomizedGreedyInitializer.o GAMatrix.o BranchAndBound.o IterativeSolver.o EMSolver.o ScaffoldExtractor.o ScaffoldComparer.o ScaffoldConverter.o GraphViz.o NWAligner.o ContigOverlapper.o optimizer.o
COBJ = Helpers.o DataStore.o DataStoreReader.o Writer.o Timers.o Reader.o ReadCoverageReader.o ReadCoverage.o ReadCoverageRepeatDetecter.o Sequence.o

include ../Makefile.config

# Need to setup prefixes
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))
_COBJ = $(patsubst %,$(CDIR)/$(ODIR)/%,$(COBJ))

.PHONY : main Common clean dirs

main : Common $(BDIR)/$(BNAME)

$(BDIR)/$(BNAME) : $(_OBJ) | $(BDIR)
	$(CCC) $(CCCFLAGS) $(CCEFLAGS) $(CCCINC) $(CCEINC) -o $(BDIR)/$(BNAME) $(_OBJ) $(_COBJ) $(CCCLIB) $(CCELIB)

dirs : | $(ODIR) $(BDIR)

$(ODIR) :
	$(MD) -p $(ODIR)

$(BDIR) :
	$(MD) -p $(BDIR)

Common :
	$(MAKE) -C $(CDIR) $(_COBJ)

$(ODIR)/%.o : %.cpp | $(ODIR)
	$(CCC) -c $(CCCFLAGS) $(CCEFLAGS) $(CCCINC) $(CCEINC) -o $@ $<

clean :
	rm -f $(BDIR)/$(BNAME)
	rm -rf $(ODIR)
