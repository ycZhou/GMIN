C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE KEYWORD
      USE commons
      USE modmxatms   ! needed for CHARMM
      USE modcharmm
      USE modamber
      USE porfuncs
      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, IX, J1, JP, NPCOUNT, NTYPEA, NPCALL, NDUMMY, INDEX, J2
      LOGICAL CAT, YESNO
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION XX, EPSAB, EPSBB, SIGAB, SIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO
      CHARACTER WORD*16
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      DOUBLE PRECISION EAMLJA0, EAMLJBETA, EAMLJZ0, DUMMY
      COMMON /EAMLJCOMM/ EAMLJA0, EAMLJBETA, EAMLJZ0
      DOUBLE PRECISION SLENGTH, EPS
      INTEGER NOK, NBAD
      COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      DOUBLE PRECISION EPS2, RAD, HEIGHT
      COMMON /CAPS/ EPS2, RAD, HEIGHT

C     LOGICAL IGNOREBIN(HISTBINMAX), FIXBIN
C     COMMON /IG/ IGNOREBIN, FIXBIN
      DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      COMMON /AMBWORD/    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      COMMON /PCALL/ NPCALL

      INTEGER NATOM, DMODE, NDUM
      DOUBLE PRECISION CHX(MXATMS), CHY(MXATMS), CHZ(MXATMS), CHMASS(MXATMS)
      CHARACTER(LEN=1) DUMMYCH
      DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN

      NPCOUNT=0
      NPCALL=0
      NSEED=0
      NS=0
      NSSTOP=0
      HIT=.FALSE.
      SAVEQ=.TRUE.
      NSAVE=5
      RESIZET=.FALSE.
      STEPOUT=.FALSE.
      SUPERSTEP=.FALSE.
      NSUPER=10
      SUPSTEP=1.1D0
      SACCRAT=0.5D0
      NSACCEPT=100
      EVSTEPT=.FALSE.
      NEVS=100
      CEIG=0.1D0
      NEVL=100
      NVECTORS=2
      TEMPS=0.8
      NRBSITES=0
      CHFREQ=1
      FTRANS=1
      FROT=1
      ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1),TEMP(1),ACCRAT(1),STEP(1),ASTEP(1),OSTEP(1),BLOCK(1),NT(1),NQ(1),EPREV(1),
     @         JUMPMOVE(1),JUMPINT(1),JDUMP(1),COORDS(3*NATOMS,1),COORDSO(3*natoms,1),VAT(natoms,1),VATO(natoms,1),
     @         JUMPTO(1))
      DO JP=1,1
         FIXSTEP(JP)=.FALSE.
         FIXTEMP(JP)=.FALSE.
         FIXBOTH(JP)=.FALSE.
         TEMP(JP)=0.3D0
         ACCRAT(JP)=0.5D0
         STEP(JP)=0.3D0
         ASTEP(JP)=0.3D0
         OSTEP(JP)=0.3D0
         BLOCK(JP)=0
         NT(JP)=0
         JUMPMOVE(JP)=.FALSE.
         JUMPINT(JP)=100
         JDUMP(JP)=.FALSE.
      ENDDO
      NEWJUMP=.FALSE.
      PNEWJUMP=0.2D0
      ECONV=0.02D0
      TABOOT=.FALSE.
      NTAB=10
      CUTOFF=1.0D6
      FINALCUTOFF=1.0D6
      MYPOWER=5
      NEON=.FALSE.
      RGCL2=.FALSE.
      AXTELL=.FALSE.
      ZSTAR=0.0D0
      GROUND=.FALSE.
      ARGON=.FALSE.
      ARNO=.FALSE.
      STAR=.FALSE.
      PLUS=.FALSE.
      TWOPLUS=.FALSE.
      DIPOLE=.FALSE.
      DUMPT=.FALSE.
      TARGET=.FALSE.
      SORTT=.FALSE.
      NTARGETS=0
      MSORIGT=.FALSE.
      MSTRANST=.FALSE.
      FRAUSIT=.FALSE.
      ANGST=.FALSE.
      MORSET=.FALSE.
      LB2T=.FALSE.
      DZTEST=.FALSE.
      ZETT1=.FALSE.
      ZETT2=.FALSE.
      P46=.FALSE.
      BLNT=.FALSE.
      DFTBT=.FALSE.
      SW=.FALSE.
      XMUL=1
      SCT=.FALSE.
      SQUEEZET=.FALSE.
      NVEC=0
C     SQUEEZER=5.0D0
C     SQUEEZED=0.95D0
      DEBUG=.FALSE.
      SEEDT=.FALSE.
      FREEZECORE=.TRUE.
      FIELDT=.FALSE.
      OHT=.FALSE.
      IHT=.FALSE.
      TDT=.FALSE.
      D5HT=.FALSE.
      CENT=.FALSE.
      FIH=0.0D0
      FTD=0.0D0
      FD5H=0.0D0
      TOLD=0.0001D0
      TOLE=0.0001D0
      CUTT=.FALSE.
      PERIODIC=.FALSE.
      NRUNS=0
      PCUT=1.0D0
      RADIUS=0.0D0
      MAXIT=500
      MAXIT2=500
      EXPFAC=10.0D0
      EXPD=1.0D0
      CQMAX=1.0D-10
      BQMAX=1.0D-3
      RHO=6.0D0
      NACCEPT=50
      NORESET=.FALSE.
      TSALLIST=.FALSE.
      QTSALLIS=0.0D0
      NPAR=1
      PARALLELT=.FALSE.
      TOSI=.FALSE.
      WELCH=.FALSE.
      BINARY=.FALSE.
      SHIFTCUT=.FALSE.
      FAL=.FALSE.
      FNI=.FALSE.
      AMBER=.FALSE.
      DPARAM=1.0D0
      FAKEWATER=.FALSE.
      AMCUT= .FALSE.
      MGBWATER=.FALSE.
      BIN=.FALSE.
      AMBERSEED= .FALSE.
      FIXT= .FALSE.
      FIX= .FALSE.
      CAP= .TRUE.
      WATERSTEP= .FALSE.
      QCUTOFF= 1.0D6
      RCUTOFF= 1.0D6
      REALQCUTOFF= 1.0D6
      REALRCUTOFF= 1.0D6
      listupdate=20

      BLJCLUSTER=.FALSE.

      CHRMMT=.FALSE.
      CHRIGIDTRANST=.FALSE.
      CHRIGIDROTT=.FALSE.
      RMST=.FALSE.
      NEWCONFT=.FALSE.
      INTMINT=.FALSE.
      DAESTAT=.FALSE.

      NOPHIPSIT=.FALSE.
      OMEGAT=.FALSE.

      BSMIN=.FALSE.
      RKMIN=.FALSE.
      PERMUTE=.FALSE.

      GAMMA=1.0D0
      TUNNELT=.FALSE.
      
      TWOD=.FALSE.
      COMPRESST=.FALSE.

      MUPDATE=4
      DGUESS=0.1D0
      BFGS=.FALSE.
      LBFGST=.TRUE.
      CONJG=.FALSE.
      TNT=.FALSE.
      TOLB=0.1D0
      DBRENTT=.FALSE.
      GUIDECUT=0.0001D0
      CPMD=.FALSE.
      DL_POLY=.FALSE.
      EFAC=0.0D0
      EAMP=0.01D0
      FIXD=.FALSE.
      NHSMOVE=1
      T12FAC=1.1D0
      RENORM=.FALSE.
      NRENORM=10
      NRENSTUCK=20
      XMOVERENORM=6.0
      TRENORM=1.0D0
      PACHECO=.FALSE.
      EAMLJT=.FALSE.
      PBGLUET=.FALSE.
      EAMALT=.FALSE.
      ALGLUET=.FALSE.
      MGGLUET=.FALSE.
      GUPTAT=.FALSE.
      WENZEL=.FALSE.
      RESTART=.FALSE.
      NEWRESTART=.FALSE.
      NRELAX=0
      NMSBSAVE=0
      AVOID=.FALSE.
      AVOIDDIST=1.0D0
      MAXSAVE=10
      NHSRESTART=0
      MAXBFGS=0.4D0

      CAPSID=.FALSE.
      STRANDT=.FALSE.
      PAHT=.FALSE.
      TIP=.FALSE.
      STOCKT=.FALSE.
      STICKYT=.FALSE.
      RIGID=.FALSE.
      TIPID=4
      HEIGHT=0.5D0
      OTPT=.FALSE.
      LJMFT=.FALSE.
      Q4T=.FALSE.
      
      THRESHOLDT=.FALSE.
      HIST=.FALSE.
      HISTRESTART=.FALSE.
      HISTSMOOTH=.FALSE.
      NSpline=1
      EPSSPHERE=0.0D0
      FIXBIN=.FALSE.

      DECAY=.FALSE.
      DECAYPARAM=0.0D0
      COOP=.FALSE.
      NCOOP=5
      COOPCUT=1.0D0

      NATBT=.FALSE.
      MAXERISE=1.0D-10
      SYMMETRIZE=.FALSE.
      NSYMINTERVAL=10
      SYMTOL1=0.1D0
      SYMTOL2=0.1D0
      SYMTOL3=0.1D0
      SYMTOL4=0.1D0
      SYMTOL5=0.1D0
      NSYMQMAX=20
      MATDIFF=0.1D0
      DISTFAC=0.0D0
      ARMA=0.4D0
      ARMB=0.4D0

      BINSTRUCTURES=.FALSE.
      TETHER=.FALSE.
      EQUIL=0
      PTMC=.FALSE.
      VISITPROP=.FALSE.
      HWINDOWS=1

      FixedEndMoveT=.FALSE.
      PIVOTP=0.0D0
      SIDECHAINP=0.0D0


      DIFFRACTT=.FALSE.
      THOMSONT=.FALSE.
      GAUSST=.FALSE.

      OPEN (5,FILE='data',STATUS='OLD')

C190   CALL INPUT(END,5)
190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF

      IF (END .OR. WORD .EQ. 'STOP') THEN

         IF (NPCOUNT.LT.NPAR) THEN
            DO J1=NPCOUNT+1,NPAR
               STEP(J1)=STEP(1)
               ASTEP(J1)=ASTEP(1)
               OSTEP(J1)=ASTEP(1)
               BLOCK(J1)=BLOCK(1)
            ENDDO
         ENDIF
        RETURN
      ENDIF

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190

      ELSE IF (WORD.EQ.'SHIFTCUT') THEN
         SHIFTCUT=.TRUE.
         CALL READF(XX)
         CUTOFF=XX

      ELSE IF (WORD.EQ.'BLJCLUSTER') THEN
         BLJCLUSTER=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
         CALL READF(CUTOFF)

      ELSE IF (WORD.EQ.'BINARY') THEN
         BINARY=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)

      ELSE IF (WORD.EQ.'SORT') THEN
         SORTT=.TRUE.

      ELSE IF (WORD.EQ.'STAR') THEN
         STAR=.TRUE.

      ELSE IF (WORD.EQ.'DIPOLES') THEN
         DIPOLE=.TRUE.

      ELSE IF (WORD.EQ.'TWOPLUS') THEN
         TWOPLUS=.TRUE.

      ELSE IF (WORD.EQ.'PLUS') THEN
         PLUS=.TRUE.

      ELSE IF (WORD.EQ.'ARGON') THEN
         ARGON=.TRUE.

      ELSE IF (WORD.EQ.'AMBER') THEN
         AMBER=.TRUE.
         CALL APARAMS
         CALL AREAD
         NATOMS=ATOMS
         DO J1=1,NATOMS
            COORDS(3*(J1-1)+1,1)=x(J1)
            COORDS(3*(J1-1)+2,1)=y(J1)
            COORDS(3*(J1-1)+3,1)=z(J1)
         ENDDO
         t=0
         ang=0
         imp=0
         count=0
      ELSE IF (WORD.EQ.'PMAX') THEN
         CALL READF(PMAX)
         WRITE(*,'(A,F14.10)') 'PMAX=  ',PMAX
      ELSE IF (WORD.EQ.'PMIN') THEN
         CALL READF(PMIN)
         WRITE(*,'(A,F14.10)') 'PMIN=  ',PMIN

      ELSE IF (WORD.EQ.'NMAX') THEN
         CALL READF(NMAX)
         WRITE(*,'(A,F14.10)') 'NMAX=  ',NMAX

      ELSE IF (WORD.EQ.'NMIN') THEN
         CALL READF(NMIN)
         WRITE(*,'(A,F14.10)') 'NMIN=  ',NMIN

      ELSE IF (WORD.EQ.'SIDESTEP') THEN
         CALL READF(SIDESTEP)
         WRITE(*,'(A,F14.10)') 'SIDESTEP=  ',SIDESTEP

      ELSE IF (WORD.EQ.'RCUTOFF') THEN
         AMCUT=.TRUE.
         CALL READF(REALRCUTOFF)
         RCUTOFF=1.1D0*REALRCUTOFF

      ELSE IF (WORD.EQ.'QCUTOFF') THEN
         AMCUT=.TRUE.
         CALL READF(REALQCUTOFF)
         QCUTOFF=1.1D0*REALQCUTOFF

      ELSE IF (WORD.EQ.'FAKEWATER') THEN
         FAKEWATER=.TRUE.
         WRITE (*,'(A)') '**********************************************************'
         WRITE (*,'(A)') '* DISTANCE DEPENDENT DIELECTRIC BEING USED - FAKE WATER! *'
         WRITE (*,'(A)') '**********************************************************'
C
C Start of CHARMM-related keywords.
C
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHRMMT=.TRUE.
C        ALLOCATE(ATMASS(NATOMS))
         CHX(1)=13.13d13 ! this way we will tell CHARMM to save it's coords into CH. arrays; othewise it will
! use input.crd only which is the default now
         CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM)
         CALL CHSETZSYMATMASS
         CALL CHALLOCATE(NATOMS)
         CALL CHSETDIHE
         IF (NATOM /= NATOMS) THEN
            WRITE(*,'(A)') 'No. of atoms in "input.crd" and file specified in CHARMM part of odata conflict'
            print *, 'NATOM,NATOMS=',NATOM, NATOMS
            call exit(10)
         ENDIF
         DO J1=1,NATOMS
            COORDS(3*(J1-1)+1,1)=CHX(J1)
            COORDS(3*(J1-1)+2,1)=CHY(J1)
            COORDS(3*(J1-1)+3,1)=CHZ(J1)
         ENDDO
         IF (INTMINT) CALL GETNINT(NINTS)  ! DJW - this is OK because CHARMM is the last keyword!
      ELSE IF (WORD.EQ.'CHARMMTYPE') THEN
         CALL READI(CHARMMTYPE)
C
C CHARMMTYPE = 1 means CHARMM22 with ACE capgroup
C CHARMMTYPE = 2 means CHARMM19 with ACE capgroup
C CHARMMTYPE = 3 means CHARMM19 no capgroup
C Used in REBUILD
C
         WRITE(*,'(A,I2)') 'CHARMMTYPE set to ',CHARMMTYPE

      ELSE IF (WORD.EQ.'CHPMAX') THEN
         CALL READF(CHPMAX)
         WRITE(*,'(A,F14.10)') 'CHPMAX=  ',CHPMAX

      ELSE IF (WORD.EQ.'CHPMIN') THEN
         CALL READF(CHPMIN)
         WRITE(*,'(A,F14.10)') 'CHPMIN=  ',CHPMIN

      ELSE IF (WORD.EQ.'CHNMAX') THEN
         CALL READF(CHNMAX)
         WRITE(*,'(A,F14.10)') 'CHNMAX=  ',CHNMAX

      ELSE IF (WORD.EQ.'CHNMIN') THEN
         CALL READF(CHNMIN)
         WRITE(*,'(A,F14.10)') 'CHNMIN=  ',CHNMIN

      ELSE IF (WORD.EQ.'NOPHIPSI') THEN
         NOPHIPSIT=.TRUE.
         WRITE(*,'(A)') 'NOPHIPSIT set : only sidechain dihedrals will be twisted'

      ELSE IF (WORD.EQ.'TOMEGA') THEN
         OMEGAT=.TRUE.
         WRITE(*,'(A)') 'TOMEGA set : peptide bonds will be twisted along with all other dihedrals'

      ELSE IF (WORD.EQ.'CHFREQ') THEN
	 CALL READI(CHFREQ)
         WRITE(*,'(A,I4,A)') 'Every ',CHFREQ,' steps TAKESTEPCH is called'

      ELSE IF (WORD.EQ.'CHRIGIDTRANS') THEN
         CHRIGIDTRANST=.TRUE.
         CALL READF(PTRANS)
         CALL READF(TRANSMAX)
         CALL READI(FTRANS)
         WRITE(*,'(A)') 'CHRIGIDTRANST set'

      ELSE IF (WORD.EQ.'CHRIGIDROT') THEN
         CHRIGIDROTT=.TRUE.
         CALL READF(PROT)
         CALL READF(ROTMAX)
         CALL READI(FROT)
         WRITE(*,'(A)') 'CHRIGIDROTT set'
  
      ELSE IF (WORD.EQ.'ARM') THEN
         ARMT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ARMA)
         IF (NITEMS.GT.2) CALL READF(ARMB)

      ELSE IF (WORD.EQ.'RMS') THEN
         RMST=.TRUE.
         CALL READF(RMSLIMIT)
         CALL READF(RMSTOL)
         CALL READI(RMSSAVE)
         CALL READI(J1)
         IF(J1.EQ.1) THEN
           SELECTT=.TRUE.
         ELSE
           SELECTT=.FALSE.
         ENDIF
         WRITE(*,'(A)') 'RMST set'
      ELSE IF (WORD.EQ.'NEWCONF') THEN
         NEWCONFT=.TRUE.
         CALL READI(NEWCONFST)
         CALL READF(NCWALL)

      ELSE IF (WORD.EQ.'INTMIN') THEN
         INTMINT=.TRUE.
C        IF (NITEMS.GT.1) THEN
C           CALL READF(IMINCUT)
C        ENDIF
C     ELSE IF (WORD.EQ.'NOCISTRANS') THEN
C        NOCISTRANS=.TRUE.

C     ELSE IF (WORD.EQ.'NORANDOM') THEN
C        NORANDOM=.TRUE.
C        IF (NITEMS.GT.1) CALL READF(RANDOMCUTOFF)

C     ELSE IF (WORD.EQ.'PERMDIHE') THEN
C        PERMDIHET=.TRUE.
C        DO J1=1,NITEMS-1
C           CALL READI(NDUM)
C           PERMDIHE(J1)=NDUM
C        ENDDO
C        NPERMDIHE=NITEMS-1
C        DO J1=1,NITEMS-1
C           print *,'PERMDIHE',PERMDIHE(J1)
C        ENDDO
C
C  End of CHARMM-related keywords.
C
      ELSE IF (WORD.EQ.'GROUND') THEN
         GROUND=.TRUE.

      ELSE IF (WORD.EQ.'NEON') THEN
         NEON=.TRUE.

      ELSE IF (WORD.EQ.'AXTELL') THEN
         AXTELL=.TRUE.
         CALL READF(ZSTAR)

      ELSE IF (WORD.EQ.'RGCL2') THEN
         RGCL2=.TRUE.

      ELSE IF (WORD.EQ.'ARNO') THEN
         ARNO=.TRUE.

      ELSE IF (WORD.EQ.'TOSI') THEN
         TOSI=.TRUE.
         CALL READF(APP)
         CALL READF(AMM)
         CALL READF(APM)
         CALL READF(RHO)

      ELSE IF (WORD.EQ.'WELCH') THEN
         WELCH=.TRUE.
         CALL READF(APP)
         CALL READF(AMM)
         CALL READF(APM)
         CALL READF(RHO)
         CALL READF(XQP)
         CALL READF(XQM)
         CALL READF(ALPHAP)
         CALL READF(ALPHAM)

      ELSE IF (WORD.EQ.'NEWJUMP') THEN
         NEWJUMP=.TRUE.
         IF (NITEMS.GT.1) CALL READF(PNEWJUMP)

      ELSE IF (WORD.EQ.'JUMPMOVE') THEN
         CALL READI(IX)
         JUMPMOVE(IX)=.TRUE.
         CALL READI(JUMPTO(IX))
         JDUMP(JUMPTO(IX))=.TRUE.
         IF (NITEMS.GT.3) CALL READI(JUMPINT(IX))

      ELSE IF (WORD.EQ.'TABOO') THEN
         TABOOT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NTAB)

      ELSE IF (WORD.EQ.'MSORIG') THEN
         MSORIGT=.TRUE.

      ELSE IF (WORD.EQ.'MSTRANS') THEN
         MSTRANST=.TRUE.

      ELSE IF (WORD.EQ.'FRAUSI') THEN
         FRAUSIT=.TRUE.

      ELSE IF (WORD.EQ.'FAL') THEN
         FAL=.TRUE.

      ELSE IF (WORD.EQ.'FNI') THEN
         FNI=.TRUE.

      ELSE IF (WORD.EQ.'ANGSTROM') THEN
         ANGST=.TRUE.

      ELSE IF (WORD.EQ.'NORESET') THEN
         NORESET=.TRUE.
C
C  Threshold acceptance rather than Metropolis, i.e. the energy change
C  can;t increase by more than a certain amount.
C
      ELSE IF (WORD.EQ.'THRESHOLD') THEN
         THRESHOLDT=.TRUE.
         PRINT*,'keyword THRESHOLD doesnt appear to do anything at the moment'
         STOP
C
C  Set Tsallis statistics with some q value.
C
      ELSE IF (WORD.EQ.'TSALLIS') THEN
         TSALLIST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(QTSALLIS)
         ENDIF

      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(ECONV)

      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
         IF (NITEMS.GT.1) CALL READF(BQMAX)

      ELSE IF (WORD.EQ.'RESIZE') THEN
         RESIZET=.TRUE.
         CALL READF(XX)
         RESIZE=XX

      ELSE IF (WORD.EQ.'DUMP') THEN
         DUMPT=.TRUE.

      ELSE IF (WORD.EQ.'TARGET') THEN
         TARGET=.TRUE.
         NTARGETS=NITEMS-1
         ALLOCATE(TARGETS(NTARGETS))
         INQUIRE(FILE='coords.target',EXIST=YESNO)
         IF (YESNO) THEN
            ALLOCATE(TCOORDS(NTARGETS,3*NATOMS))
            OPEN(UNIT=1,FILE='coords.target',STATUS='OLD')
            READ(1,*) ((TCOORDS(J1,J2),J2=1,3*NATOMS),J1=1,NTARGETS)
            CLOSE(1)
         ENDIF
         DO J1=2,NITEMS
            CALL READF(XX)
            TARGETS(J1-1)=XX
         ENDDO

      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)
      ELSE IF (WORD.EQ.'MAXIT') THEN
         CALL READI(IX)
         MAXIT=IX
         IF (NITEMS.GT.2) THEN
            CALL READI(IX)
            MAXIT2=IX
         ENDIF

      ELSE IF (WORD.EQ.'ACCEPTRATIO') THEN
         IF (NITEMS-1.GT.NPAR) THEN
            PRINT '(A)','Number of acceptance ratios exceeds NPAR - quit'
            STOP
         ENDIF
         DO J1=1,NITEMS-1
            CALL READF(ACCRAT(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               ACCRAT(J1)=ACCRAT(1)
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(IX)
         NACCEPT=IX

      ELSE IF (WORD.EQ.'PERIODIC') THEN
         PERIODIC=.TRUE.
         CALL READF(XX)
         BOXLX=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            BOXLY=XX
            CALL READF(XX)
            BOXLZ=XX
         ELSE
            BOXLY=BOXLX
            BOXLZ=BOXLX
         ENDIF

      ELSE IF (WORD.EQ.'STEPS') THEN
         NRUNS=1
         CALL READI(IX)
         MCSTEPS(1)=IX
         CALL READF(XX)
         TFAC(1)=XX
         IF (NITEMS.GT.3) THEN
            NRUNS=2
            CALL READI(IX)
            MCSTEPS(2)=IX
            CALL READF(XX)
            TFAC(2)=XX
         ENDIF
         IF (NITEMS.GT.5) THEN
            NRUNS=3
            CALL READI(IX)
            MCSTEPS(3)=IX
            CALL READF(XX)
            TFAC(3)=XX
         ENDIF

      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.

      ELSE IF (WORD.EQ.'RADIUS') THEN
         CALL READF(XX)
         RADIUS=XX

      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF

      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
         CALL READF(CQMAX)
      
      ELSE IF (WORD.EQ.'FIXBOTH') THEN
         IF (NITEMS.EQ.1) THEN
            FIXBOTH(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXBOTH(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FIXTEMP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXTEMP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXTEMP(IX)=.TRUE.
            ENDDO
         ENDIF
      ELSE IF (WORD.EQ.'FIXSTEP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXSTEP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXSTEP(IX)=.TRUE.
            ENDDO
         ENDIF
C
C Read in the maximum initial step size, factor for determining angular
C moves, and for rigid bodies the angular step size and the size of the
C blocks for Cartesian and angular moves.
C
C For parallel runs different values can be used for different runs by
C adding additional "STEP" lines to the data file. Otherwise the
C parameters for subsequent parallel runs are set to the values for the
C first one.
C
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            PRINT '(A)','Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))
         IF (NITEMS.GT.3) CALL READF(OSTEP(NPCOUNT))
         IF (NITEMS.GT.4) CALL READI(BLOCK(NPCOUNT))

      ELSE IF (WORD.EQ.'CUTOFF') THEN
         CUTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(CUTOFF)
         FINALCUTOFF=CUTOFF
         IF (NITEMS.GT.2) CALL READF(FINALCUTOFF)

      ELSE IF (WORD.EQ.'SEED') THEN
         SEEDT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(IX)
            NSSTOP=IX
         ENDIF

      ELSE IF (WORD.EQ.'NOFREEZE') THEN
         FREEZECORE=.FALSE.

      ELSE IF (WORD.EQ.'P46') THEN
         P46=.TRUE.
         BLNT=.TRUE.

      ELSE IF (WORD.EQ.'BLN') THEN
         BLNT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         READ(1,*) DUMMYCH
         READ(1,*) LJREPBB, LJATTBB
         READ(1,*) LJREPLL, LJATTLL
         READ(1,*) LJREPNN, LJATTNN
         READ(1,*) DUMMYCH
         READ(1,*) DUMMYCH
         READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(1)
         PRINT '(A,I8,A)','BLN sequence of ',NATOMS,' beads read:'
         WRITE(*,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         PRINT '(A)',' '
         PRINT '(A,I8,A)','BLN dihedral types:'
         WRITE(*,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         PRINT '(A)',' '
         PRINT '(A,2F15.5)','B-B LJ coefficients: ',LJREPBB, LJATTBB
         PRINT '(A,2F15.5)','L-L LJ coefficients: ',LJREPLL, LJATTLL
         PRINT '(A,2F15.5)','N-N LJ coefficients: ',LJREPNN, LJATTNN
         PRINT '(A,4F15.5)','Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         PRINT '(A,4F15.5)','Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         PRINT '(A,4F15.5)','Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
C        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
C    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS) 

      ELSE IF (WORD.EQ.'DFTB') THEN
         DFTBT=.TRUE.

      ELSE IF (WORD.EQ.'SW') THEN
         SW=.TRUE.

      ELSE IF (WORD.EQ.'MULTIPLICITY') THEN
         CALL READI(XMUL)

      ELSE IF (WORD.EQ.'MORSE') THEN
         MORSET=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(XX)
            RHO=XX
         ENDIF

      ELSE IF (WORD.EQ.'LB2') THEN
         LB2T=.TRUE.

      ELSE IF (WORD.EQ.'DZUGUTOV') THEN
         DZTEST=.TRUE.
         CALL READF(DZP1)
         CALL READF(DZP2)
         CALL READF(DZP3)
         CALL READF(DZP4)
         CALL READF(DZP5)
         CALL READF(DZP6)
         CALL READF(DZP7)

      ELSE IF (WORD.EQ.'ZETT1') THEN
         ZETT1=.TRUE.

      ELSE IF (WORD.EQ.'ZETT2') THEN
         ZETT2=.TRUE.

      ELSE IF (WORD.EQ.'SC') THEN
         SCT=.TRUE.
         CALL READI(IX)
         NN=IX
         CALL READI(IX)
         MM=IX
         CALL READF(XX)
         SIG=XX
         CALL READF(XX)
         SCEPS=XX
         CALL READF(XX)
         SCC=XX

      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
      ELSE IF (WORD.EQ.'PACHECO') THEN
         PACHECO=.TRUE.
      ELSE IF (WORD.EQ.'EAMAL') THEN
         EAMALT=.TRUE.
      ELSE IF (WORD.EQ.'PBGLUE') THEN
         PBGLUET=.TRUE.
      ELSE IF (WORD.EQ.'ALGLUE') THEN
         ALGLUET=.TRUE.
      ELSE IF (WORD.EQ.'MGGLUE') THEN
         MGGLUET=.TRUE.
      ELSE IF (WORD.EQ.'GUPTA') THEN
         GUPTAT=.TRUE.
         CALL READI(GATOM)
      ELSE IF (WORD.EQ.'EAMLJ') THEN
         EAMLJT=.TRUE.
         CALL READF(EAMLJA0)
         CALL READF(EAMLJBETA)
         CALL READF(EAMLJZ0)
      ELSE IF (WORD.EQ.'COMPRESS') THEN
         COMPRESST=.TRUE.
         CALL READF(COMP)
C
C  MYPOWER provides a means to set the initial premultiplication factor for the
C  gradient in MYLINMIN
C
      ELSE IF (WORD.EQ.'POWER') THEN
         CALL READI(IX)
         MYPOWER=IX
C
C PARALLEL must come before STEP and ACCRAT
C
      ELSE IF (WORD.EQ.'PARALLEL') THEN
         PARALLELT=.TRUE.
         CALL READI(NPAR)
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,JUMPMOVE,JUMPINT,JDUMP,COORDS,NQ,
     @              JUMPTO,EPREV,COORDSO,VAT,VATO) 
         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR),OSTEP(NPAR), 
     @         BLOCK(NPAR),NT(NPAR),JUMPMOVE(NPAR),JUMPINT(NPAR),JDUMP(NPAR),NQ(NPAR),JUMPTO(NPAR),COORDS(3*NATOMS,NPAR),
     @         COORDSO(3*natoms,NPAR),VAT(natoms,NPAR),VATO(natoms,NPAR),EPREV(NPAR))
         DO JP=1,NPAR
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
            OSTEP(JP)=0.3D0
            BLOCK(JP)=0
            NT(JP)=0
            JUMPMOVE(JP)=.FALSE.
            JUMPINT(JP)=100
            JDUMP(JP)=.FALSE.
         ENDDO
      ELSE IF (WORD.EQ.'2D') THEN
         TWOD=.TRUE.

C
C  Number of BFGS updates before resetting, default=4
C
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
C
C  Initial guess for diagonal matrix elements in LBFGS.
C
       ELSE IF (WORD.EQ.'DGUESS') THEN
          CALL READF(DGUESS)
      ELSE IF (WORD.EQ.'BFGS') THEN
         BFGS=.TRUE.
      ELSE IF (WORD.EQ.'TN') THEN
         TNT=.TRUE.
         PRINT '(A)','optimisation with tn no longer supported'
         STOP
      ELSE IF (WORD.EQ.'DBRENT') THEN
         DBRENTT=.TRUE.
      ELSE IF (WORD.EQ.'TOLBRENT') THEN
         CALL READF(TOLB)
C
C  Conjugate gradient optimisation instead of default LBFGS
C
      ELSE IF (WORD.EQ.'CG') THEN
         LBFGST=.FALSE.
         CONJG=.TRUE.
      ELSE IF (WORD.EQ.'DIELEC') THEN
         CALL READF(XX)
         DPARAM=XX
         WRITE(*,'(A,F9.5)') ' Dielectric constant = ',DPARAM
      ELSE IF (WORD.EQ.'GUIDE') THEN
         CALL READF(GUIDECUT)
      ELSE IF (WORD.EQ.'STICKY') THEN
         STICKYT=.TRUE.
         RIGID=.TRUE.
         CALL READI(NRBSITES)
         CALL READF(STICKYSIG)
         PRINT*,'NRBSITES=',NRBSITES 
         PRINT*,'STICKYSIG=',STICKYSIG 
         ALLOCATE(SITE(NRBSITES,3))
         DO J1=1,NRBSITES
            READ(5,*) SITE(J1,1:3)
C           CALL READF(SITE(J1,1))
C           CALL READF(SITE(J1,2))
C           CALL READF(SITE(J1,3))
            PRINT '(A,I5,3G20.10)','J1,site: ',J1,SITE(J1,1:3)
         ENDDO
      ELSE IF (WORD.EQ.'TIP') THEN
         TIP=.TRUE.
         RIGID=.TRUE.
         IF (NITEMS.GT.1) CALL READI(TIPID)
         IF (TIPID.EQ.5) NRBSITES=5
         IF (TIPID.EQ.4) NRBSITES=4
         IF (TIPID.EQ.3) NRBSITES=3
         IF (TIPID.EQ.2) NRBSITES=4
         IF (TIPID.EQ.1) NRBSITES=3
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD.EQ.'STOCK') THEN
         STOCKT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(STOCKMU)
         CALL READF(STOCKLAMBDA)
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD.EQ.'STRAND') THEN
         STRANDT=.TRUE.
         RIGID=.TRUE.
C
C  The nine reference site positions per strand.
C
         NRBSITES=9
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=-2.7298862082
         SITE(1,2)=2.3622865625 
         SITE(1,3)=0.6475151629
         SITE(2,1)=-1.7492122114
         SITE(2,2)=2.3331194664 
         SITE(2,3)=0.5887015133
         SITE(3,1)=-1.5963638586
         SITE(3,2)=1.4304320585 
         SITE(3,3)=0.2442792479
         SITE(4,1)=-0.6166461313
         SITE(4,2)=1.4301805389 
         SITE(4,3)=0.1327546571
         SITE(5,1)=-0.4460267836
         SITE(5,2)=0.5254809645  
         SITE(5,3)=-0.2196837962
         SITE(6,1)=0.5313983749 
         SITE(6,2)=0.5210707739  
         SITE(6,3)=-0.3409645197
         SITE(7,1)=0.7065341613  
         SITE(7,2)=-0.3914277962 
         SITE(7,3)=-0.6719579835 
         SITE(8,1)=1.6776397940  
         SITE(8,2)=-0.3830053500 
         SITE(8,3)=-0.8355266604
         SITE(9,1)=1.8162689403  
         SITE(9,2)=-1.3093381947 
         SITE(9,3)=-1.1427874015
      ELSE IF (WORD.EQ.'PAH') THEN
         PAHT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=36
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD.EQ.'CAPSID') THEN
         CAPSID=.TRUE.
         RIGID=.TRUE.
         CALL READF(RHO)
         CALL READF(EPS2)
         CALL READF(RAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)
C
C  The six reference site positions per capped pentagon. These need to
C  be multiplied by RAD, including the repulsive site!
C
         NRBSITES=6
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=1.0D0
         SITE(1,2)=0.0D0
         SITE(1,3)=0.0D0
         SITE(2,1)=((-1.0D0 + Sqrt(5.0D0)))/4.0D0
         SITE(2,2)=(Sqrt((5.0D0 + Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(2,3)=0.0D0
         SITE(3,1)=((-1.0D0 - Sqrt(5.0D0)))/4.0D0
         SITE(3,2)=(Sqrt((5.0D0 - Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(3,3)=0.0D0
         SITE(4,1)=((-1 - Sqrt(5.0D0)))/4.0D0
         SITE(4,2)=-(Sqrt((5.0D0 - Sqrt(5.0D0))/2.))/2.0D0
         SITE(4,3)=0.0D0
         SITE(5,1)=((-1 + Sqrt(5.0D0)))/4.0D0
         SITE(5,2)=-(Sqrt((5.0D0 + Sqrt(5.0D0))/2.))/2.0D0
         SITE(5,3)=0.0D0
         SITE(6,1)=0.0D0
         SITE(6,2)=0.0D0
         SITE(6,3)=HEIGHT
      ELSE IF (WORD.EQ.'CPMD') THEN
         CPMD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') ' ERROR - no CPMD system specified'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 12
            ENDIF
         ENDDO
12       CONTINUE
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)
      ELSE IF (WORD.EQ.'BSMIN') THEN
         BSMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
      ELSE IF (WORD.EQ.'RKMIN') THEN
         RKMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
      ELSE IF (WORD.EQ.'PERMUTE') THEN
         PERMUTE=.TRUE.
C
C  Undocumented options (some probably obsolete/broken).
C
      ELSE IF (WORD.EQ.'EVSTEP') THEN
         EVSTEPT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NEVS)
         IF (NITEMS.GT.2) CALL READF(CEIG)
         IF (NITEMS.GT.3) CALL READI(NEVL)
         IF (NITEMS.GT.4) CALL READI(NVECTORS)

      ELSE IF (WORD.EQ.'SUPERSTEP') THEN
         SUPERSTEP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSUPER)
         IF (NITEMS.GT.2) CALL READF(SUPSTEP)
         IF (NITEMS.GT.3) CALL READF(TEMPS)
         IF (NITEMS.GT.4) CALL READF(SACCRAT)
         IF (NITEMS.GT.5) CALL READI(NSACCEPT)

      ELSE IF (WORD.EQ.'STEPOUT') THEN
         STEPOUT=.TRUE.

      ELSE IF (WORD.EQ.'LJMF') THEN
         LJMFT=.TRUE.
C        CALL LJPARAMMF
         PRINT '(A)','LJMF not currently maintained'
         STOP
      ELSE IF (WORD.EQ.'OTP') THEN
         OTPT=.TRUE.
         RIGID=.TRUE.
C        CALL OTPPARAMMF
      ELSE IF (WORD.EQ.'WENZEL') THEN
         WENZEL=.TRUE.

C
C  Take hard sphere type moves.
C  T12FAC is the fraction of the first collision time to be used in HSMOVE
C
      ELSE IF (WORD.EQ.'FIXD') THEN
         FIXD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NHSMOVE)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(T12FAC)
         ENDIF
C
C  Renormalisation attempt
C
C  TRENORM is the temperature for the Metropolis accept/reject comparison
C          of lowest energies calculated over NRENORM steps having moved
C          XMOVERENORM atoms. NRENORM is dynamically adjusted with a
C          minimum value equal to half the original NRENORM.
C
      ELSE IF (WORD.EQ.'RENORM') THEN
         RENORM=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRENORM)
         IF (NITEMS.GT.2) CALL READF(XMOVERENORM)
         IF (NITEMS.GT.3) CALL READF(TRENORM)
         IF (NITEMS.GT.4) CALL READI(NRENSTUCK)

      ELSE IF (WORD.EQ.'GAMMA') THEN
         CALL READF(GAMMA)
         TUNNELT=.TRUE.

      ELSE IF (WORD .EQ. 'TOLD') THEN
        CALL READF(XX)
        TOLD=XX

      ELSE IF (WORD .EQ. 'TOLE') THEN
        CALL READF(XX)
        TOLE=XX

C     ELSE IF (WORD.EQ.'SQUEEZE') THEN
C        CALL READI(NVEC)
C        SQUEEZET=.TRUE.
C        IF (NITEMS.GT.2) THEN
C           CALL READF(XX)
C           SQUEEZER=XX
C        ENDIF
C        IF (NITEMS.GT.3) THEN
C           CALL READF(XX)
C           SQUEEZED=XX
C        ENDIF
      ELSE IF (WORD.EQ.'OH') THEN
         FIELDT=.TRUE.
         OHT=.TRUE.
         CALL READF(XX)
         FOH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(EXPFAC)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF

      ELSE IF (WORD.EQ.'IH') THEN
         FIELDT=.TRUE.
         IHT=.TRUE.
         CALL READF(XX)
         FIH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'TD') THEN
         FIELDT=.TRUE.
         TDT=.TRUE.
         CALL READF(XX)
         FTD=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'D5H') THEN
         FIELDT=.TRUE.
         D5HT=.TRUE.
         CALL READF(XX)
         FD5H=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF
      ELSE IF (WORD.EQ.'DL_POLY') THEN
         DL_POLY=.TRUE.
         CALL READI(NATOMS)
C
C  Reseed runs if a step is not accepted in twice the relaxation time,
C  defined in terms of a number of mc steps NRELAX. NHSRESTART defines
C  the number of hard sphere moves used to produce the new starting
C  configuration. 
C
      ELSE IF (WORD.EQ.'RESTART') THEN
         RESTART=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
C
C  Reseed runs if the energy does not decrease within NRELAX mc steps.
C  NHSRESTART defines the number of hard sphere moves used to produce the new starting
C  configuration. If NHSRESTART=0 then the geometry is changed using RESEED.
C
      ELSE IF (WORD.EQ.'NEWRESTART') THEN
         NEWRESTART=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
         IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
         IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
C
C  Read data for previous geometries that lead to reseeding, which
C  probably approximate MSB bottoms.
C
      ELSE IF (WORD.EQ.'READMSB') THEN
         INQUIRE(FILE='MSBdata',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(*,'(A)') 'ERROR - READMSB specified, but no MSBdata data file found'
            STOP
         ELSE
            IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
            IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
            OPEN(UNIT=34,FILE='MSBdata',STATUS='OLD')
57          READ(34,*,END=56) DUMMY
            NMSBSAVE=NMSBSAVE+1
            MSBE(NMSBSAVE)=DUMMY
            READ(34,*) (MSBCOORDS(J1,NMSBSAVE),J1=1,3*NATOMS)
            IF (NMSBSAVE.LT.MAXSAVE) GOTO 57
56          WRITE(*,'(A,I6,A)') 
     1         'Energies and coordinates read for ',NMSBSAVE,' previous structures from MSBdata'
            CLOSE(34)
         ENDIF
C
C  Specify resetting if the latest structure gets too close to minima saved
C  in MSBCOORDS. Use bipartite matching and closest approach distance AVOIDDIST.
C  Maximum number of saved strutures is specified by MAXSAVE.
C 
      ELSE IF (WORD.EQ.'AVOID') THEN
         AVOID=.TRUE.
         IF (NITEMS.GT.1) CALL READF(AVOIDDIST)
         IF (NITEMS.GT.2) CALL READI(MAXSAVE)
         IF (.NOT.ALLOCATED(MSBCOORDS)) THEN
            ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
         ELSE
            PRINT*,'reallocating MSBCOORDS'
            DEALLOCATE(MSBCOORDS)
            ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
         ENDIF
         IF (.NOT.ALLOCATED(MSBE)) THEN
            ALLOCATE(MSBE(MAXSAVE))
         ELSE
            PRINT*,'reallocating MSBE'
            DEALLOCATE(MSBE)
            ALLOCATE(MSBE(MAXSAVE))
         ENDIF
      ELSE IF (WORD.EQ.'EXPFAC') THEN
         CALL READF(EFAC)
         IF (NITEMS.GT.2) CALL READF(EAMP)

      ELSE IF (WORD.EQ.'FIXEDEND') THEN
         FixedEndMoveT = .TRUE.
         IF (NITEMS.GT.1) CALL READF(PIVOTP)
         IF (NITEMS.GT.2) CALL READF(SIDECHAINP)


C tvb -=== Basin-sampling related keywords: ===-

      ELSE IF (WORD.EQ.'HISTOGRAM') THEN
         HIST=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTINT)
         CALL READF(HISTFAC)
         CALL READI(HBINS)
         CALL READF(HISTFACMUL)
         CALL READI(TargetWL)
         CALL READF(HPERCENT)
         ALLOCATE(HDIST(HBINS),HWEIGHT(HBINS),HISTVALS(HBINS),LHISTVALS(HBINS),IGNOREBIN(HBINS))
         DO J1=1,HBINS
            HISTVALS(J1)=0
            LHISTVALS(J1)=0
            HWEIGHT(J1)=1.0D0
            HDIST(J1)=0.0D0
C           DO J2=1,HBINS
C              HTRANS(J1,J2)=1.0D0 ! transition matrix
C           ENDDO
         ENDDO
C
C  During the run HDIST contains the sum of the distances found for minima in each bin. The
C  average is saved in hist.new.
C
            DO J1=1,HBINS
               HDIST(J1)=HDIST(J1)*HISTVALS(J1)
            ENDDO

      ELSE IF (WORD.EQ.'HISTRESTART') THEN
         HISTRESTART=.TRUE.

      ELSE IF (WORD.EQ.'HISTSMOOTH') THEN
         CALL READI(NSpline)
C Parameters of the temperature range on which to calculate thermodynamic properties in Basin Sampling

      ELSE IF (WORD.EQ.'HISTTEMP') THEN
         CALL READF(MinimalTemperature)
         CALL READF(MaximalTemperature)
         CALL READI(NTempPoints)

C Saves every n'th structure to the file with corresponding bin label
      ELSE IF (WORD.EQ.'BINSTRUCTURES') THEN
         BINSTRUCTURES=.TRUE.
         CALL READI(SaveNth)
C Tethered WL walk to determine anharmonic vibrational density of states
      ELSE IF (WORD.EQ.'TETHER') THEN
         TETHER=.TRUE.
         CALL READF(hdistconstraint)
         CALL READI(hwindows)
         lhbins=int(hbins/hwindows)
         CALL READF(ExtrapolationPercent)
         lhbins=int(hbins/hwindows)
         sampledbins=int((1.0d0-ExtrapolationPercent)*hbins/hwindows)
         CALL READF(lnHarmFreq)


C Accumulation of thermodynamic statistics starting after Equil steps, calculated thermodynamic properties is dumped every DumpEveryNthQuench quench.
      ELSE IF (WORD.EQ.'EQUILIBRATION') THEN
         CALL READI(EQUIL)
         CALL READI(DumpEveryNthQuench)

C Choice of convergence regime: histogram flatness (default), 
C VisitProp - minimal number of visits proportional to 1/sqrt(ln(f))
      ELSE IF (WORD.EQ.'VISITPROP') THEN
         VISITPROP=.TRUE.
C
C Parallel tempering Monte Carlo
C 'PT parameters ', NumberOfTrajectories, MinimalTemperature, MaximalTemperature, MCCycles, MCCyclesEquilibration 
      ELSE IF (WORD.EQ.'PTMC') THEN
         PTMC=.TRUE.
         CALL READI(NumberOfTrajectories)
         CALL READF(MinimalTemperature)
         CALL READF(MaximalTemperature)
         CALL READI(MCCycles)
         CALL READI(MCCyclesEquilibration)
         CALL READF(SwapProb)

C Request calculation of structural order parameter Q4 on the fly 
      ELSE IF (WORD.EQ.'Q4') THEN
         Q4T=.TRUE.

C tvb -=== End of basin-sampling related keywords: ===-

C
C  Correlated random moves, in the sense that the magnitude of the step
C  decays exponentially as DECAYPARAM from a randomly chosen atom.
C
      ELSE IF (WORD.EQ.'DECAY') THEN
         DECAY=.TRUE.
         CALL READF(DECAYPARAM)


C
C  Alternative correlated moves: NCOOP nearest neighbours of a randomly selected atom
C  all move by the same amount. NCOOP default is 5.
C
      ELSE IF (WORD.EQ.'COOPMOVE') THEN
         COOP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCOOP)
         IF (NITEMS.GT.2) CALL READF(COOPCUT)
      ELSE IF (WORD.EQ.'NATB') THEN
         NATBT=.TRUE.
      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)
C
C  Keyword and parameters for symmetrisation.
C
      ELSE IF (WORD.EQ.'SYMMETRISE') THEN
         SYMMETRIZE=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSYMINTERVAL)
         IF (NITEMS.GT.2) CALL READF(SYMTOL1)
         IF (NITEMS.GT.3) CALL READF(SYMTOL2)
         IF (NITEMS.GT.4) CALL READF(SYMTOL3)
         IF (NITEMS.GT.5) CALL READF(SYMTOL4)
         IF (NITEMS.GT.6) CALL READF(SYMTOL5)
         IF (NITEMS.GT.7) CALL READI(NSYMQMAX)
         IF (NITEMS.GT.8) CALL READF(MATDIFF) ! appears to have little effect now
         IF (NITEMS.GT.9) CALL READF(DISTFAC)
      ELSE IF (WORD.EQ.'THOMSON') THEN
         THOMSONT=.TRUE.
      ELSE IF (WORD.EQ.'GAUSS') THEN
         GAUSST=.TRUE.
         CALL READI(GMODES) ! number of nodes
      ELSE IF (WORD.EQ.'DIFFRACT') THEN
         DIFFRACTT=.TRUE.
      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF
      CALL FLUSH(6)
      GOTO 190

      RETURN
      END
