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
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/FFTs.o \
	${OBJECTDIR}/cluster_analysis.o \
	${OBJECTDIR}/fitting_planes.o \
	${OBJECTDIR}/functions.o \
	${OBJECTDIR}/global_variables.o \
	${OBJECTDIR}/integration_methods.o \
	${OBJECTDIR}/interpolation_fitting_methods.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/nrtype.o \
	${OBJECTDIR}/nrutil.o \
	${OBJECTDIR}/rheology_routines.o \
	${OBJECTDIR}/sinal_processing.o \
	${OBJECTDIR}/statistical_functions.o \
	${OBJECTDIR}/swarm_routines.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analise_estru_dinamica

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analise_estru_dinamica: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analise_estru_dinamica ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/FFTs.o: FFTs.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/FFTs.o FFTs.F90

${OBJECTDIR}/cluster_analysis.o: cluster_analysis.f90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/cluster_analysis.o cluster_analysis.f90

${OBJECTDIR}/fitting_planes.o: fitting_planes.f90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/fitting_planes.o fitting_planes.f90

${OBJECTDIR}/functions.o: functions.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/functions.o functions.F90

${OBJECTDIR}/global_variables.o: global_variables.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/global_variables.o global_variables.F90

${OBJECTDIR}/integration_methods.o: integration_methods.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/integration_methods.o integration_methods.F90

${OBJECTDIR}/interpolation_fitting_methods.o: interpolation_fitting_methods.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/interpolation_fitting_methods.o interpolation_fitting_methods.F90

${OBJECTDIR}/main.o: main.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/main.o main.F90

${OBJECTDIR}/nrtype.o: nrtype.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/nrtype.o nrtype.F90

${OBJECTDIR}/nrutil.o: nrutil.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/nrutil.o nrutil.F90

${OBJECTDIR}/rheology_routines.o: rheology_routines.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/rheology_routines.o rheology_routines.F90

${OBJECTDIR}/sinal_processing.o: sinal_processing.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/sinal_processing.o sinal_processing.F90

${OBJECTDIR}/statistical_functions.o: statistical_functions.f90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/statistical_functions.o statistical_functions.f90

${OBJECTDIR}/swarm_routines.o: swarm_routines.F90
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/swarm_routines.o swarm_routines.F90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
