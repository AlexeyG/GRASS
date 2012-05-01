CCC = g++

CCCFLAGS = -O2 -m64 -std=gnu++0x -Wall -fPIC -fexceptions -fopenmp -DNDEBUG -DIL_STD
CCCINC = -I../Common/ -I/data/bio/alexeygritsenk/apps/include -I/data/bio/alexeygritsenk/apps/ILOG/cplex/include -I/data/bio/alexeygritsenk/apps/ILOG/concert/include -I/data/bio/alexeygritsenk/apps/include/ncbi-tools++
CCCLIB = -L/data/bio/alexeygritsenk/apps/lib -lbamtools -lxalgoalignnw -lxobjmgr -lgenome_collection -lseqset -lseqedit -lseq -lseqcode -lsequtil -lpub -lmedline -lbiblio -lgeneral -lxser -lxutil -lxncbi -ltables -L/data/bio/alexeygritsenk/apps/ILOG/cplex/lib/x86-64_sles10_4.1/static_pic -lilocplex -lcplex -L/data/bio/alexeygritsenk/apps/ILOG/concert/lib/x86-64_sles10_4.1/static_pic -lconcert -lm -pthread

MAKE = make
MD = mkdir

CDIR = ../Common
ODIR = obj
BDIR = ../bin

BNAME = scaffoldOptimizer
OBJ = Configuration.o OverlapperConfiguration.o DPGraph.o DPSolver.o MIQPSolver.o GAIndividual.o GASolver.o FixedMIQPSolver.o ExtendedFixedMIQPSolver.o RelaxedFixedMIQPSolver.o SolverConfiguration.o RandomizedGreedyInitializer.o GAMatrix.o BranchAndBound.o IterativeSolver.o EMSolver.o ScaffoldExtractor.o ScaffoldComparer.o ScaffoldConverter.o GraphViz.o NWAligner.o ContigOverlapper.o optimizer.o
COBJ = Helpers.o DataStore.o DataStoreReader.o Writer.o Timers.o Reader.o ReadCoverageReader.o ReadCoverage.o ReadCoverageRepeatDetecter.o Sequence.o

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
