include Makefile.config

all : main extra

main : Common breakpointCounter coverageUtil dataLinker scaffoldOptimizer

extra : dataFilter dataLinker dataSelector dataSimulator kmer readCleaner readDiff

Common :
	$(MAKE) -C Common -f $(CDIR)/Makefile

breakpointCounter :
	$(MAKE) -C breakpointCounter -f $(CDIR)/Makefile

coverageUtil :
	$(MAKE) -C coverageUtil -f $(CDIR)/Makefile

dataLinker :
	$(MAKE) -C dataLinker -f $(CDIR)/Makefile

scaffoldOptimizer :
	$(MAKE) -C breakpointCounter -f $(CDIR)/Makefile

clean :
	$(RM) -rf $(ODIR)
