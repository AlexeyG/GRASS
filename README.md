=GRASS: a generic algorithm for scaffolding next-generation sequencing assemblies=

This SVN repository contains source code implementing algorithms used in "GRASS: a generic algorithm for scaffolding next-generation sequencing assemblies", Gritsenko et al. (in preparation).

The follwing modules are available:
  * Scaffold optimizer (solves the MIQP optimization problem formulation to produce scaffold nucleotide sequences)
    * Data linker (uses available information sources to derive scaffolding constraints; currently supports paired ends and mate pairs, and reference genomes)

    Additionally the following tools are present:
      * Breakpoint counter (assesses scaffold correctness by aligning scaffolds to reference sequences and counting breakpoints)
