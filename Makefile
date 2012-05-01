include Makefile.config

all : main extra

main : Common breakpointCounter coverageUtil dataLinker scaffoldOptimizer

extra : dataFilter dataSelector dataSimulator kmer readCleaner readDiff

Common :
	$(MAKE) -C Common

breakpointCounter :
	$(MAKE) -C breakpointCounter core

coverageUtil :
	$(MAKE) -C coverageUtil core

dataLinker :
	$(MAKE) -C dataLinker core

scaffoldOptimizer :
	$(MAKE) -C scaffoldOptimizer core

dataFilter :
	$(MAKE) -C dataFilter core

dataSelector :
	$(MAKE) -C dataSelector core

dataSimulator :
	$(MAKE) -C dataSimulator core

kmer :
	$(MAKE) -C kmer core

readCleaner :
	$(MAKE) -C readCleaner core

readDiff :
	$(MAKE) -C readDiff core

clean :
	$(RM) -rf bin
	$(RM) -rf */obj
