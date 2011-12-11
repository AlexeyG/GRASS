#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/GAIndividual.o \
	${OBJECTDIR}/ExtendedFixedMIQPSolver.o \
	${OBJECTDIR}/EMSolver.o \
	${OBJECTDIR}/DPGraph.o \
	${OBJECTDIR}/SolverConfiguration.o \
	${OBJECTDIR}/MIQPSolver.o \
	${OBJECTDIR}/FixedMIQPSolver.o \
	${OBJECTDIR}/Configuration.o \
	${OBJECTDIR}/RelaxedFixedMIQPSolver.o \
	${OBJECTDIR}/GAMatrix.o \
	${OBJECTDIR}/BranchAndBound.o \
	${OBJECTDIR}/optimizer.o \
	${OBJECTDIR}/ScaffoldExtractor.o \
	${OBJECTDIR}/_ext/740233504/NWAligner.o \
	${OBJECTDIR}/ScaffoldComparer.o \
	${OBJECTDIR}/DPSolver.o \
	${OBJECTDIR}/GASolver.o \
	${OBJECTDIR}/GraphViz.o \
	${OBJECTDIR}/RandomizedGreedyInitializer.o \
	${OBJECTDIR}/IterativeSolver.o \
	${OBJECTDIR}/ScaffoldConverter.o \
	${OBJECTDIR}/Solver.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/scaffoldoptimizer

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/scaffoldoptimizer: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/scaffoldoptimizer ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/GAIndividual.o: GAIndividual.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GAIndividual.o GAIndividual.cpp

${OBJECTDIR}/ExtendedFixedMIQPSolver.o: ExtendedFixedMIQPSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ExtendedFixedMIQPSolver.o ExtendedFixedMIQPSolver.cpp

${OBJECTDIR}/EMSolver.o: EMSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/EMSolver.o EMSolver.cpp

${OBJECTDIR}/DPGraph.o: DPGraph.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/DPGraph.o DPGraph.cpp

${OBJECTDIR}/SolverConfiguration.o: SolverConfiguration.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/SolverConfiguration.o SolverConfiguration.cpp

${OBJECTDIR}/MIQPSolver.o: MIQPSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/MIQPSolver.o MIQPSolver.cpp

${OBJECTDIR}/FixedMIQPSolver.o: FixedMIQPSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/FixedMIQPSolver.o FixedMIQPSolver.cpp

${OBJECTDIR}/Configuration.o: Configuration.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/Configuration.o Configuration.cpp

${OBJECTDIR}/RelaxedFixedMIQPSolver.o: RelaxedFixedMIQPSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/RelaxedFixedMIQPSolver.o RelaxedFixedMIQPSolver.cpp

${OBJECTDIR}/GAMatrix.o: GAMatrix.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GAMatrix.o GAMatrix.cpp

${OBJECTDIR}/BranchAndBound.o: BranchAndBound.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/BranchAndBound.o BranchAndBound.cpp

${OBJECTDIR}/optimizer.o: optimizer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/optimizer.o optimizer.cpp

${OBJECTDIR}/ScaffoldExtractor.o: ScaffoldExtractor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ScaffoldExtractor.o ScaffoldExtractor.cpp

${OBJECTDIR}/_ext/740233504/NWAligner.o: /Users/alexeyg/Documents/src/tud-scaffolding/scaffoldOptimizer/NWAligner.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/740233504
	${RM} $@.d
	$(COMPILE.c) -g -I../Common -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/740233504/NWAligner.o /Users/alexeyg/Documents/src/tud-scaffolding/scaffoldOptimizer/NWAligner.cpp

${OBJECTDIR}/ScaffoldComparer.o: ScaffoldComparer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ScaffoldComparer.o ScaffoldComparer.cpp

${OBJECTDIR}/DPSolver.o: DPSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/DPSolver.o DPSolver.cpp

${OBJECTDIR}/GASolver.o: GASolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GASolver.o GASolver.cpp

${OBJECTDIR}/GraphViz.o: GraphViz.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/GraphViz.o GraphViz.cpp

${OBJECTDIR}/RandomizedGreedyInitializer.o: RandomizedGreedyInitializer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/RandomizedGreedyInitializer.o RandomizedGreedyInitializer.cpp

${OBJECTDIR}/IterativeSolver.o: IterativeSolver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/IterativeSolver.o IterativeSolver.cpp

${OBJECTDIR}/ScaffoldConverter.o: ScaffoldConverter.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/ScaffoldConverter.o ScaffoldConverter.cpp

${OBJECTDIR}/Solver.o: Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Common -I/Users/alexeyg/apps/include -I/Users/alexeyg/apps/include/ncbi-tools++ -I/Users/alexeyg/apps/ILOG/cplex/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/Solver.o Solver.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/scaffoldoptimizer

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
