GPROG =	GMIN.2.2

OBJS1 = commons.o countatoms.o modamb.o modcharmm.o modmxatms.o
OBJS2 =	porfuncs.o centre.o finalio.o modconsts_trans_97.o modconsts.o dist.o\
	io1.o keyword.o main.o mc.o mcruns.o morse.o \
	potential.o quench.o rad.o dprand.o saveit.o seed.o \
	sort.o sort2.o sort3.o sort4.o takestep.o mycpu_time.o trans.o \
	finalq.o symmetry.o ptgrp.o eigsrt.o SiSW.o taboo.o reseed.o newinertia.o supermc.o \
	tosifumi.o ortho.o compress.o mylbfgs.o input.o ddfpmin.o dlnsrch.o cgmin.o linmin.o \
	brent.o mnbrak.o dbrent.o f1dim.o zwischen.o hsmove.o PachecoC60.o AT.o EAMLJ_sub.o \
	mind.o Pbglue.o wenzel.o odesd.o capsid.o rigidfuncs.o tip.o pah.o strand.o \
	SW.o qmod.o ljpbin.o fdm.o dftb.o ljpshift.o dzugutov.o \
        fd.o fedor.o welch.o BLJcluster.o stock.o Farkas.o \
        sc.o Zetterling.o MSorig.o MSorigc.o MStrans.97.o \
	frausi.o p46mer.o lj.o modperm.o modf1com.o mododesd.o amber.o EAMal.o Alglue.o Mgglue.o \
	Gupta.o orient.o Natb.o sticky.o enumerate.o minperm.o minpermdist.o LB2.o \
	dgetrf.o dgetri.o reorient.o thomson.o Q4.o basinsampling.o tether.o tetherfunc.o BLN.o diffract.o \
        gauss.o
OBJS = ${OBJS1} ${OBJS2}
GENF90FILES = porfuncs.f90
CHDUM = chdummy.o
VPATH = .:CHARMM
LDFLAGS = -L.


###################################### COMPILERS AND COMPILER FLAGS ######################################
#
###### The Portland Group Compiler Technology Fortran 90 compiler
# FC = pgf90
# FFLAGS= -Mextend -fastsse -O3 -Mvect=sse,assoc,cachesize:1024000,recog,transform 
# FFLAGS= -Mextend -O3 -Munroll -Mscalarsse -Mnoframe -Mvect=sse -Mcache_align -Mflushz -tp k8-32 -Msecond_underscore
# FFLAGS= -Mextend -C -g
#
#  NOOPT = -O0
#  SEARCH_PATH =  -I..
# the double reference to -lblas seems to be needed!
#  LIBS = -lblas -llapack -lblas 
#  LDFLAGS= -L.
##### end of The Portland Group Compiler Technology Fortran 90 compiler
#
##### NAGWare Fortran 95 compiler flags (SAT)
#FC = f95
#FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O4
# FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte
# NOOPT = -O0
#SEARCH_PATH = -I..
#LDFLAGS= -L.
#SWITCH=nag
# the double reference to -lblas seems to be needed!
#LIBS = -lblas -llapack -lblas
##### end of NAGWare Fortran 95 compiler flags
#
##########################################################
#
##### Intel compilers
# 
# FC = ifort
# FC = ifc
# FFLAGS= -132 -C -g
# FFLAGS= -132 -O3 
# FFLAGS= -132 -Vaxlib -O3 # for ifc
# FFLAGS= -132 -O3 -ip -static # -ipo or -fast doesnt work 
# NOOPT= -O0
# SWITCH=ifort
# SEARCH_PATH =  -I..
# LDFLAGS= -L.
# LIBS = -lblas -llapack -lblas -lblas -llapack
##### end of Intel compiler flags
#
##### Sun Studio compilers
#
# FC = f90
# FFLAGS = -e -C -g 
# SEARCH_PATH = -M..
# FFLAGS = -e -native -libmil -fsimple=2 -dalign -xlibmopt -depend -fnonstd -xO5 \
#          -xtarget=ultra3cu -xcache=64/32/4:8192/512/2 -xprefetch=yes -xvector=yes -xcrossfile # -XlistE -Xlistwar -z muldefs
# LIBS = -xlic_lib=sunperf -L/usr/SUNWspro/lib
# LDFLAGS= -L.
# NOOPT=
# the double reference to -lblas seems to be needed!
# LIBS = -lblas -llapack -lblas 

###################################### RULES AND TARGETS ######################################
.SUFFIXES:
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<

default: $(GPROG)

CGMIN1: $(CPROG1)

CGMIN2: $(CPROG2)

$(GPROG): $(OBJS) $(EXTRAS) $(CHDUM)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS) $(OBJS) $(CHDUM) $(LDFLAGS) $(LIBS)


clean:
	rm -f ${GMIN} ${CGMIN1} ${CGMIN2} $(OBJS) *.mod $(EXTRAS) *.lst *.o pref.dat prefx.msg ifc?????? nag.f90 FOR021.DAT porfuncs.f90
	cd BLAS; make clean
	cd LAPACK; make clean

timing:
	rm -f GMIN.2.0
	$(FC) $(FFLAGS) -p $(OBJS) $(EXTRAS) -o $(GPROG) $(LIBS)

feedback:
	$(FC) $(FFLAGS) -xprofile=use:gmin.profile *.f -o $(GPROG)

rebuild:
	make clean
	make

libblas.a: SAT-Ghost
	cd BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}"
liblapack.a: SAT-Ghost
	cd LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}"
SAT-Ghost:

porfuncs.f90: porfuncs.csh
	./porfuncs.csh ${SWITCH} > porfuncs.f90

###################################### DEPENDENCIES ######################################
${GPROG}: libblas.a liblapack.a
${CGPROG1}: libblas.a liblapack.a libcharmm.a
${CGPROG2}: libblas.a liblapack.a libcharmm.a
libcharmm.a: commons.o modcharmm.o modmxatms.o

${OBJS2}: ${OBJS1}

Alglue.o:      commons.o
BLJcluster.o:  commons.o
EAMLJ_sub.o:   commons.o
EAMal.o:       commons.o
Gupta.o:       commons.o
MSorig.o:      commons.o modconsts.o dist.o
MSorigc.o:     commons.o modconsts.o dist.o
MStrans.97.o:  commons.o modconsts_trans_97.o dist.o
Mgglue.o:      commons.o
PachecoC60.o:  commons.o
Pbglue.o:      commons.o
SW.o:          commons.o
SiSW.o:        commons.o
Zetterling.o:  commons.o
amber.o:       commons.o modamb.o
bmin.o:        commons.o modamb.o
capsid.o:      commons.o
centre.o:      commons.o
cgmin.o:       commons.o
compress.o:    commons.o
compress2.o:   commons.o
ddfpmin.o:     commons.o
dftb.o:        commons.o
dzugutov.o:    commons.o
evstep.o:      commons.o
f1dim.o:       commons.o modf1com.o
fd.o:          commons.o
fdm.o:         commons.o
finalio.o:     commons.o modamb.o qmod.o modcharmm.o
finalq.o:      commons.o qmod.o
frausi.o:      dist.o
grnd.o:        commons.o
hmat1n.o:      commons.o
hmatd_.o:      commons.o
hsmove.o:      commons.o
io1.o:         commons.o modamb.o modperm.o qmod.o modcharmm.o porfuncs.o
keyword.o:     commons.o modamb.o modcharmm.o porfuncs.o
linmin.o:      commons.o modf1com.o
lj.o:          commons.o
ljpbin.o:      commons.o
ljpshift.o:    commons.o
main.o:        commons.o modf1com.o countatoms.o modperm.o qmod.o modamb.o
mc.o:          qmod.o modcharmm.o porfuncs.o
mcruns.o:      commons.o
mind.o:        porfuncs.o
mf.o:          commons.o
mnbrak.o:      commons.o
commons.o:  countatoms.o
morse.o:       commons.o
mylbfgs.o:     commons.o modamb.o porfuncs.o
odesd.o:       commons.o mododesd.o
olami.o:       commons.o
otp.o:         commons.o
p46mer.o:      commons.o
BLN.o:      commons.o
pah.o:         commons.o
potential.o:   commons.o modperm.o qmod.o modcharmm.o porfuncs.o
ptgrp.o:       commons.o
quench.o:      commons.o qmod.o porfuncs.o
rad.o:         commons.o
rdpot.o:       commons.o
saveit.o:      commons.o qmod.o
sc.o:          commons.o
seed.o:        commons.o
strand.o:      commons.o
supermc.o:     commons.o
symmetry.o:    commons.o porfuncs.o
taboo.o:       commons.o
takestep.o:    commons.o
tip.o:         commons.o
tosifumi.o:    commons.o
welch.o:       commons.o
zwischen.o:    commons.o modf1com.o 
stock.o:       commons.o
sticky.o:      commons.o
diffract.o:    NP.H
diffract.o:    commons.o
tether.o:      tetherfunc.o
mycpu_time.o:  commons.o
gauss.o:       commons.o
