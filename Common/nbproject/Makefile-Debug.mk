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
	${OBJECTDIR}/_ext/152718998/DataStoreReader.o \
	${OBJECTDIR}/_ext/152718998/ReadCoverageWriter.o \
	${OBJECTDIR}/_ext/152718998/DataStore.o \
	${OBJECTDIR}/_ext/152718998/ReadCoverageReader.o \
	${OBJECTDIR}/_ext/152718998/Writer.o \
	${OBJECTDIR}/_ext/152718998/Helpers.o \
	${OBJECTDIR}/_ext/152718998/ReadCoverageRepeatDetecter.o \
	${OBJECTDIR}/_ext/152718998/Aligner.o \
	${OBJECTDIR}/_ext/152718998/Reader.o \
	${OBJECTDIR}/_ext/152718998/XATag.o \
	${OBJECTDIR}/_ext/152718998/AlignmentReader.o \
	${OBJECTDIR}/_ext/152718998/AlignerConfiguration.o \
	${OBJECTDIR}/_ext/152718998/ReadCoverage.o \
	${OBJECTDIR}/_ext/152718998/Converter.o \
	${OBJECTDIR}/_ext/152718998/Sequence.o \
	${OBJECTDIR}/_ext/152718998/DataStoreWriter.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCommon.dylib

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCommon.dylib: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -dynamiclib -install_name libCommon.dylib -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCommon.dylib -fPIC ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/152718998/DataStoreReader.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStoreReader.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/DataStoreReader.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStoreReader.cpp

${OBJECTDIR}/_ext/152718998/ReadCoverageWriter.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageWriter.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/ReadCoverageWriter.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageWriter.cpp

${OBJECTDIR}/_ext/152718998/DataStore.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStore.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/DataStore.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStore.cpp

${OBJECTDIR}/_ext/152718998/ReadCoverageReader.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageReader.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/ReadCoverageReader.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageReader.cpp

${OBJECTDIR}/_ext/152718998/Writer.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Writer.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Writer.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Writer.cpp

${OBJECTDIR}/_ext/152718998/Helpers.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Helpers.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Helpers.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Helpers.cpp

${OBJECTDIR}/_ext/152718998/ReadCoverageRepeatDetecter.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageRepeatDetecter.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/ReadCoverageRepeatDetecter.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverageRepeatDetecter.cpp

${OBJECTDIR}/_ext/152718998/Aligner.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Aligner.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Aligner.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Aligner.cpp

${OBJECTDIR}/_ext/152718998/Reader.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Reader.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Reader.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Reader.cpp

${OBJECTDIR}/_ext/152718998/XATag.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/XATag.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/XATag.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/XATag.cpp

${OBJECTDIR}/_ext/152718998/AlignmentReader.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/AlignmentReader.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/AlignmentReader.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/AlignmentReader.cpp

${OBJECTDIR}/_ext/152718998/AlignerConfiguration.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/AlignerConfiguration.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/AlignerConfiguration.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/AlignerConfiguration.cpp

${OBJECTDIR}/_ext/152718998/ReadCoverage.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverage.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/ReadCoverage.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/ReadCoverage.cpp

${OBJECTDIR}/_ext/152718998/Converter.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Converter.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Converter.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Converter.cpp

${OBJECTDIR}/_ext/152718998/Sequence.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/Sequence.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/Sequence.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/Sequence.cpp

${OBJECTDIR}/_ext/152718998/DataStoreWriter.o: /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStoreWriter.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/152718998
	${RM} $@.d
	$(COMPILE.cc) -g -I/Users/alexeyg/apps/include -fPIC  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/152718998/DataStoreWriter.o /Users/alexeyg/Documents/src/tud-scaffolding/Common/DataStoreWriter.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCommon.dylib

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
