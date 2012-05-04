include Makefile.config

.PHONY : main extra all Common breakpointCounter coverageUtil dataLinker scaffoldOptimizer dataFilter dataSelector dataSimulator kmer readCleaner readDiff tar manual

main : Common breakpointCounter scaffoldOptimizer 

extra : dataFilter dataSelector dataSimulator kmer readCleaner readDiff coverageUtil

all : main extra

Common :
	$(MAKE) -C Common

breakpointCounter : Common
	$(MAKE) -C breakpointCounter

coverageUtil : Common
	$(MAKE) -C coverageUtil

dataLinker : Common
	$(MAKE) -C dataLinker

scaffoldOptimizer : Common
	$(MAKE) -C scaffoldOptimizer

dataFilter : Common
	$(MAKE) -C dataFilter

dataSelector : Common
	$(MAKE) -C dataSelector

dataSimulator : Common
	$(MAKE) -C dataSimulator

kmer : Common
	$(MAKE) -C kmer

readCleaner : Common
	$(MAKE) -C readCleaner

readDiff : Common
	$(MAKE) -C readDiff

manual :
	$(MAKE) -C manual

clean :
	$(RM) -rf bin
	$(RM) -rf */obj
	$(MAKE) -C manual clean

tar :
	tar -cf grass-src.tar LICENSE Makefile Makefile.config */*.cpp */*.h */Makefile manual/manual.tex

