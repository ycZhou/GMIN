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
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
      USE commons
      USE modcharmm
      USE porfuncs
      IMPLICIT NONE
      

      INTEGER J1, NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), J2, NSTEPS, JP, 
     1        UNT, ITERATIONS, NSUPERCOUNT, NQTOT, JACCPREV, NREN, NLAST, NSTEPREN, BRUN,QDONE,JBEST(NPAR),
     2        NRMS
      INTEGER :: NSYMCALL=0
      DOUBLE PRECISION POTEL, SCALEFAC, RANDOM, DPRAND, 
     1                 TIME, SPOTEL(NSUPER), SCOORDS(3*NATOMS,NSUPER), SCREENC(3*NATOMS),
     2                 EPPREV(NPAR), QSTART, QFINISH, RANNJ, RMIN, RMINO, RCOORDS(3*NATOMS),ELASTSYM(NPAR),
     3                 RCOORDSO(3*NATOMS), RVAT(NATOMS), RVATO(NATOMS), EPSSAVE, EBEST(NPAR),
     4                 BESTCOORDS(3*NATOMS,NPAR), endtime, RMSD
      LOGICAL CHANGEDE
      CHARACTER FNAME*9

      LOGICAL EVAP, ATEST, STAY, evapreject
      COMMON /EV/ EVAP, evapreject
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      COMMON /Q4C/ QSTART, QFINISH

      ALLOCATE(TMOVE(NPAR), OMOVE(NPAR))
      NSTEPREN=0
      evapreject=.FALSE.

!      IF (PTMC) THEN
!         CALL ParallelTempering
!         RETURN
!      ENDIF


C tvb requesting a basin-sampling MC run:
      
      IF (HIST.and.(.not.TETHER)) then
         CALL BasinSampling
         RETURN
      ELSEIF (TETHER) THEN
         CALL TetheredWL
         RETURN
      ENDIF
C

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      NRMS=0
      NLAST=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM
      DO JP=1,NPAR
         TMOVE(JP)=.TRUE.
         OMOVE(JP)=.TRUE.
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
         IF (JDUMP(JP).AND.(.NOT.NEWJUMP)) THEN
            WRITE(FNAME,'(A,I1)') 'ebuffer.',JP
            UNT=70+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
            WRITE(FNAME,'(A,I1)') 'cbuffer.',JP
            UNT=70+NPAR+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
         ENDIF
      ENDDO
C
C  Calculate the initial energy and save in EPREV
C
      WRITE(*,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
         IF (NPAR.GT.1) THEN
            WRITE(*,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,'] Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ELSE
            WRITE(*,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ENDIF
C
C  EPREV saves the previous energy in the Markov chain.
C  EBEST and JBEST record the lowest energy since the last reseeding and the
C  step it was attained at. BESTCOORDS contains the corresponding coordinates.
C
         EPREV(JP)=POTEL
         EPPREV(JP)=0.0D0
         ELASTSYM(JP)=0.0D0
         EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=0
         RMIN=POTEL
         RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
         CALL FLUSH(6)
         COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
      ENDDO
      EPSSPHERE=EPSSAVE

      IF (NPAR.EQ.1) THEN
         WRITE(*,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      ELSE
         WRITE(*,'(A,I3,A,I10,A)') 'Starting ',NPAR,' parallel MC runs of ',NSTEPS,' steps'
      ENDIF
      WRITE(*,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

      NSUPERCOUNT=NSUPER
C
C  *********************************** Main basin-hopping loop *********************************
C
      DO J1=1,NSTEPS
         IF (NEWJUMP) RANNJ=DPRAND()
C
C  ********************************* Loop over NPAR parallel runs ******************************
C
         DO JP=1,NPAR 
            IF (RIGID.AND.(BLOCK(JP).GT.0)) THEN
               IF (MOD(J1-1,BLOCK(JP)).EQ.0) THEN
                  IF (MOD((J1-1)/BLOCK(JP),2).EQ.0) THEN
                     WRITE(*,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body translational moves'
                     TMOVE(JP)=.TRUE.
                     OMOVE(JP)=.FALSE.
                  ELSE 
                     WRITE(*,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body angular moves'
                     OMOVE(JP)=.TRUE.
                     TMOVE(JP)=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
C
C  Jump-moves.
C
            IF (JUMPMOVE(JP).AND.(MOD(J1,JUMPINT(JP)).EQ.0)) CALL JUMPM(RANNJ,J1,JP,EPPREV)
C
C  Ordinary steps.
C
23          CONTINUE

! Don;t call symmetry unless the minimum in the Markov chain has changed.
! We should really check if the minimum has changed since the last call to SYMMETRY,
! which can be done with ABS(ELASTSYM(JP)-EPREV(JP)) if NSYMINTERVAL=1.

C           IF (SYMMETRIZE.AND.(MOD(J1-1,NSYMINTERVAL).EQ.0).AND.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) THEN
            IF ((.NOT.SEEDT).AND.SYMMETRIZE.AND.(MOD(J1-1,NSYMINTERVAL).EQ.0)) THEN
               IF ((ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) NSYMCALL=0
C              PRINT '(A,3G20.10,I5)','ELASTSYM,EPREV,diff,NSYMCALL=',ELASTSYM(JP),EPREV(JP),ABS(ELASTSYM(JP)-EPREV(JP)),NSYMCALL
               ELASTSYM(JP)=EPREV(JP)
               CALL SYMMETRY(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL)
               IF (HIT) GOTO 37
C
C  Check for reseeding.
C
               POTEL=EPREV(JP) ! NEWRES assumes POTEL is the energy of the current structure in COORDS
!              PRINT '(A,3F15.5,L10)','in mc 1 EPREV(JP),POTEL,EBEST,CHANGEDE=',EPREV(JP),POTEL,EBEST,CHANGEDE
               IF (CHANGEDE.AND.NEWRESTART) THEN
                  CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV)
               ENDIF
!              PRINT '(A,3F15.5)','in mc 2 EPREV(JP),POTEL,EBEST=',EPREV(JP),POTEL,EBEST
            ELSEIF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
               NSYMREM=0                                       ! Should therefore reset NSYMREM.
            ENDIF

!           ELSE
               IF (AMBER) THEN
                  CALL TAKESTEPAM(JP)
               ELSE IF (CHRMMT) THEN
                 IF(CHRIGIDTRANST.AND.MOD(J1,FTRANS).EQ.0) 
     &             CALL MKRIGIDTRANS(JP)
                 IF(CHRIGIDROTT.AND.MOD(J1,FROT).EQ.0) 
     &             CALL MKRIGIDROT(JP)
                 IF(MOD(J1,CHFREQ).EQ.0) CALL TAKESTEPCH(JP)
               ELSE
                  CALL TAKESTEP(JP)
               ENDIF
               NQ(JP)=NQ(JP)+1
               CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
               NQTOT=NQTOT+1
C
C  Check for results of taboo list. SELFT is set in taboo.
C
               IF (SELFT) THEN
                  CALL RESET(JP,NATOMS,NPAR,NSEED,COORDS,COORDSO,VAT,VATO)
                  IF (DEBUG) THEN
                     PRINT*,'Taboo list:'
                     WRITE(*,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
                     WRITE(*,73) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
73                   FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' TABOO')
                  ENDIF
                  GOTO 23
               ENDIF
   
!              PRINT '(A,3F15.5)','in mc 3 EPREV(JP),POTEL,EBEST=',EPREV(JP),POTEL,EBEST
               IF (NPAR.GT.1) THEN
                  WRITE(*,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,'] Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
               ELSE
                  WRITE(*,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
               ENDIF
!           ENDIF
            CALL FLUSH(6)
C
C  RMS compared to reference structure 'compare'
            IF(RMST.AND.CHRMMT) THEN
              CALL CHRMS(JP,RMSD)
              WRITE(*,'(A,F15.5)')'RMSD = ',RMSD
              IF(RMSD.LE.RMSLIMIT) CALL SAVERMS(JP,POTEL,RMSD)
C                NRMS=NRMS+1
C                WRITE(CNRMS,'(I6)') NRMS
C                CALL CHARMMDUMP(COORDS,'rms.'//TRIM(ADJUSTL(CNRMS)))
C                OPEN(UNIT=20,FILE='rms.'//TRIM(ADJUSTL(CNRMS)),POSITION='APPEND',STATUS='OLD')
C                WRITE(20,'(A,I6,A,F15.5,A,F15.5)') '*   Qu ',NQ(JP),' E=',POTEL,' RMSD=',RMSD
C                CLOSE(20)
C              ENDIF
            ENDIF
C
C DAESTAT keyword prints two files :
C stat.all which contains all quenches and their energies
C stat.acc which contains only accepted quenches and their energies
C used to analyse how many new minima are being found
C
            IF (DAESTAT) THEN
              PRINT*,'DAESTAT block in mc.f not implemented'
              STOP
C              CALL CALCMIND(JP,MIND)
C              CALL CALCDIHE(DIHE)
C              WRITE(*,'(A,I6,3F20.10)') 'NQALL POTEL',NQ(JP),POTEL,MIND,DIHE
C              WRITE(36,'(I6,3F20.10)') NQ(JP),POTEL,MIND,DIHE
            ENDIF
C
C  Check for reseeding.
C
            IF (NEWRESTART.AND.(.NOT.SEEDT)) CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV)
            IF (STAY) THEN
            ELSE IF (EVAP .and. .not.evapreject) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL RESET(JP,NATOMS,NPAR,NSEED,COORDS,COORDSO,VAT,VATO)
               IF (DEBUG) THEN
                  WRITE(*,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE 
C
C Accept or reject step. If the quench did not converge then allow a
C potenial move, but count it as a rejection in terms of NSUCCESS and
C NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
C However, for TIP5P some cold fusion events that had not actually reached the threshold for
C rejection were actually accepted. Must prevent this!
C
               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
               ELSE
                  CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM)
               ENDIF
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(*,34) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV(JP)).GT.ECONV) THEN
C                    NRELAX=J1-JACCPREV
C                    IF (RESTART) WRITE(*,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
                     JACCPREV=J1
                  ENDIF
                  IF (QDONE.EQ.1) THEN
                     NSUCCESS(JP)=NSUCCESS(JP)+1
                  ELSE
                     NFAIL(JP)=NFAIL(JP)+1
                  ENDIF
                  EPPREV(JP)=EPREV(JP)
                  EPREV(JP)=POTEL
                  COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                  VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               ELSE
                  NFAIL(JP)=NFAIL(JP)+1
                  CALL RESET(JP,NATOMS,NPAR,NSEED,COORDS,COORDSO,VAT,VATO)
                  IF (DEBUG) THEN
                     WRITE(*,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
C           WRITE(*,'(A,4F20.10)') 'Q4 values ',QS,QF,QSTART,QFINISH
C
C  Dump coordinates and energy if another run is attempting to jump to this one.
C
            IF ((NPAR.GT.1).AND.(.NOT.NEWJUMP)) CALL DUMPJ(JP,JUMPTO,NPAR,COORDS,NATOMS,EPREV)
C
C  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
C
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
C
C  Check the acceptance ratio.
C
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)

            TEMP(JP)=TEMP(JP)*SCALEFAC
            IF (HIT) GOTO 37
         ENDDO
C
C  ****************************** End of loop over NPAR parallel runs *****************************
C
C  Supersteps
C
         IF (SUPERSTEP) CALL SUPERMC(SPOTEL,SCOORDS,NSUPERCOUNT,POTEL)
C
C  Renormalised type steps
C
         IF (RENORM) CALL REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
         IF (STAY) CALL RESET(1,NATOMS,NPAR,NSEED,COORDS,COORDSO,VAT,VATO)
  
         CALL FLUSH(6)
      ENDDO
C
C  ******************************* End of main basin-hopping loop *******************************
C

37    CONTINUE
      DO JP=1,NPAR
         WRITE(*,10) JP,NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
10       FORMAT('[',I2,'] Acceptance ratio for run=',F12.5,' Step=',F12.5,
     1          ' Angular step factor=',F12.5,' Temperature=',F12.5)
         IF (RIGID) WRITE(*,11)
11       FORMAT('[',I2,'] Rigid body orientational step=',F12.5)
      ENDDO
     
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM)
      USE commons
      use qmodule
      IMPLICIT NONE
      DOUBLE PRECISION ENEW, EOLD, DUMMY, DPRAND, RANDOM, EREF, TEOLD, TENEW, RATIO
      DOUBLE PRECISION TRANS, DISTMIN, DISTMINOLD
      LOGICAL ATEST, FLAT, evap, evapreject
      INTEGER NP,INDEXOLD, INDEXNEW, J1, NDUMMY
      DATA DISTMINOLD /0.0D0/
      COMMON /DMIN/ DISTMIN
C     COMMON /IG/ IGNOREBIN, FIXBIN
      common /ev/ evap, evapreject

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
      IF (TUNNELT) THEN
         TEOLD=TRANS(EOLD,QMIN(1),GAMMA)
         TENEW=TRANS(ENEW,QMIN(1),GAMMA)
C        WRITE(*,'(A,4F20.10)') 'TEOLD,TENEW,QMIN(1),GAMMA=',TEOLD,TENEW,QMIN(1),GAMMA
      ELSE
         TEOLD=EOLD
         TENEW=ENEW
      ENDIF

      IF (TSALLIST) THEN
         EREF=QMIN(NP)*1.1D0
         DUMMY=(1.0D0-(1.0D0-QTSALLIS)*(TENEW-EREF)/TEMP(NP))/(1.0D0-(1.0D0-QTSALLIS)*(TEOLD-EREF)/TEMP(NP))
         DUMMY=DUMMY**(QTSALLIS/(1.0D0-QTSALLIS))
C        WRITE(*,'(A,4F20.10)') 'TENEW,TEOLD,EREF,DUMMY=',TENEW,TEOLD,EREF,DUMMY
         IF (DUMMY.GE.1.0D0) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DUMMY.GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF
      ELSE
C
C  Standard canonical sampling.
C
         IF (TENEW.LT.TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(TENEW-TEOLD)/MAX(TEMP(NP),1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF
      ENDIF 

      RETURN 
      END 

      SUBROUTINE JUMPM(RANNJ,J1,JP,EPPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, J2, UNT, NDUM, ITERATIONS, BRUN, QDONE
      DOUBLE PRECISION RANNJ, RANDOM, DPRAND, EPPREV(NPAR), DUMMY, TIME, EJ, SCREENC(3*NATOMS)

      IF (NEWJUMP) THEN
         IF (RANNJ.LT.PNEWJUMP) THEN
            RANDOM=DPRAND()
            IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
C           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))/TEMP(JP)).GT.RANDOM) THEN
C           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
               WRITE(*,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' accepted before quench ',J1
               DUMMY=EPREV(JP)
               EPREV(JP)=EPREV(JUMPTO(JP))
               EPREV(JUMPTO(JP))=DUMMY
               DUMMY=EPPREV(JP)
               EPPREV(JP)=EPPREV(JUMPTO(JP))
               EPPREV(JUMPTO(JP))=DUMMY
               DO J2=1,NATOMS
                  DUMMY=VATO(J2,JP)
                  VATO(J2,JP)=VATO(J2,JUMPTO(JP))
                  VATO(J2,JUMPTO(JP))=DUMMY
                  DUMMY=VAT(J2,JP)
                  VAT(J2,JP)=VAT(J2,JUMPTO(JP))
                  VAT(J2,JUMPTO(JP))=DUMMY
               ENDDO
               DO J2=1,3*NATOMS
                  DUMMY=COORDS(J2,JP)
                  COORDS(J2,JP)=COORDS(J2,JUMPTO(JP))
                  COORDS(J2,JUMPTO(JP))=DUMMY
                  DUMMY=COORDSO(J2,JP)
                  COORDSO(J2,JP)=COORDSO(J2,JUMPTO(JP))
                  COORDSO(J2,JUMPTO(JP))=DUMMY
               ENDDO
            ELSE
               WRITE(*,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                  ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' rejected before quench ',J1
            ENDIF
         ENDIF
      ELSE
         UNT=70+JUMPTO(JP)
         REWIND(UNT)           
         NDUM=INT(DPRAND()*(NQ(JUMPTO(JP))-1))
         IF (DEBUG) PRINT*,'Should be choosing buffer energy number ',NDUM
         DO J2=1,NDUM
            READ(UNT,*) EJ
         ENDDO
         DO J2=1,NQ(JUMPTO(JP))-1-NDUM
            READ(UNT,*) DUMMY
         ENDDO
         RANDOM=DPRAND()
C
C  Coordinates are only read if the jump is successful.
C
         IF (DEXP((EPREV(JP)-EJ)*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
            WRITE(*,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' accepted before quench ',NQ(JP)
            EPREV(JP)=EJ
            UNT=70+NPAR+JUMPTO(JP)
            REWIND(UNT)           
            DO J2=1,(NDUM-1)*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
            READ(UNT,*) (COORDS(J2,JP),J2=1,3*NATOMS)
C
C  Coordinates should be converged already, but we need to reset VAT and VATO.
C
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            DO J2=1,NATOMS
               VATO(J2,JP)=VAT(J2,JP)
            ENDDO
            DO J2=1,3*NATOMS
               COORDSO(J2,JP)=COORDS(J2,JP)
            ENDDO
            IF (DEBUG) THEN
               PRINT*,'Jump coordinates:'
               WRITE(*,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
            ENDIF
            DO J2=1,NATOMS*(NQ(JUMPTO(JP))-1)-NDUM*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
         ELSE
            WRITE(*,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' rejected before quench ',NQ(JP)
         ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE DUMPJ(JP,JUMPTO,NPAR,COORDS,NATOMS,EPREV)
      IMPLICIT NONE
      LOGICAL TEST
      INTEGER NPAR, J2, JP, JUMPTO(NPAR), NATOMS, UNT
      DOUBLE PRECISION COORDS(3*NATOMS,NPAR), EPREV(NPAR)

      DO J2=1,NPAR
         IF (JUMPTO(J2).EQ.JP) TEST=.TRUE.
      ENDDO
C     TEST=.FALSE.
C     IF (TEST) THEN
C        UNT=70+JP
C        WRITE(UNT,'(F20.10)') EPREV(JP)
C        UNT=70+NPAR+JP
C        WRITE(UNT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
C     ENDIF

      RETURN
      END

      SUBROUTINE RESET(JP,NATOMS,NPAR,NSEED,COORDS,COORDSO,VAT,VATO)
      IMPLICIT NONE
      INTEGER JP, NSEED, J2, NATOMS, NPAR
      DOUBLE PRECISION COORDS(3*NATOMS,NPAR), COORDSO(3*NATOMS,NPAR), VAT(NATOMS,NPAR), VATO(NATOMS,NPAR)

      DO J2=1,3*(NATOMS-NSEED)
         COORDS(J2,JP)=COORDSO(J2,JP)
      ENDDO
      DO J2=1,NATOMS
         VAT(J2,JP)=VATO(J2,JP)
      ENDDO
  
      RETURN
      END
C
C
C
      SUBROUTINE REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
      USE commons
      IMPLICIT NONE
      INTEGER ITERATIONS, J2, JACCPREV, J1, NQTOT, BRUN, QDONE
      DOUBLE PRECISION TIME, POTEL, RCOORDS(3*NATOMS), RMIN, RVAT(NATOMS), SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

10    CALL HSMOVE(COORDS,1,NHSRESTART)
      CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
C
C  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
C
      IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
         DO J2=1,3*NATOMS
            COORDS(J2,1)=COORDSO(J2,1)
         ENDDO
         GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(*,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(*,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ(1),' E=',
     1              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
      DO J2=1,3*NATOMS
         COORDSO(J2,1)=COORDS(J2,1)
         RCOORDS(J2)=COORDS(J2,1)
      ENDDO
      DO J2=1,NATOMS
         VATO(J2,1)=VAT(J2,1)
         RVAT(J2)=VAT(J2,1)
      ENDDO
      EPREV(1)=POTEL
      RMIN=POTEL

      RETURN
      END
C
C
C
      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      USE commons
      USE modcharmm
      IMPLICIT NONE
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), JP, NFAILT(NPAR), NSUCCESST(NPAR), J1, J2, NDUMMY
      LOGICAL evap, evapreject
      DOUBLE PRECISION DUMMY, DUMMY2, DUMMY3, DUMMY4, HWMAX,P0,FAC
C     COMMON /IG/ IGNOREBIN, FIXBIN
C     COMMON /MOVE/ TMOVE, OMOVE
      common /ev/ evap, evapreject

      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
      
      IF (P0.GT.ACCRAT(JP)) THEN
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.05D0
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)/1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=NHSMOVE+1 
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)*1.05D0
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)*1.05D0
               ELSE
                  STEP(JP)=FAC*STEP(JP)
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=FAC*TRANSMAX
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=FAC*ROTMAX  
               ENDIF
            ENDIF
            ASTEP(JP)=ASTEP(JP)*1.05D0
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)/1.05D0
         ENDIF
      ELSE
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.D0/1.05D0
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)*1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=MAX(1,NHSMOVE-1)
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)/1.05D0
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/1.05D0
               ELSE
                  STEP(JP)=FAC*STEP(JP)
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=FAC*TRANSMAX
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=FAC*ROTMAX
               ENDIF
            ENDIF
            ASTEP(JP)=ASTEP(JP)/1.05D0
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)*1.05D0
         ENDIF
      ENDIF
C
      WRITE(*,'(A,I2,A,I4,A,F8.4,A,F8.4)')
     1  '[',JP,'] Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      IF (FIXBOTH(JP)) THEN
      ELSE IF (FIXSTEP(JP)) THEN
         IF(.NOT.FIXTEMP(JP)) WRITE(*,'(A,I2,A,F12.4)') '[',JP,'] Temperature is now:',TEMP(JP)
      ELSE
         WRITE(*,'(A,I2,A)',ADVANCE='NO') '[',JP,'] Steps are now:'
         WRITE(*,'(A,F8.4)',ADVANCE='NO') '  STEP=',STEP(JP)    
         IF(ASTEP(JP).GT.0.D0) WRITE(*,'(A,F8.4)',ADVANCE='NO')'  ASTEP=',ASTEP(JP) 
         IF(CHRIGIDTRANST.AND.CHRMMT) WRITE(*,'(A,F8.4)',ADVANCE='NO')'  TRANSMAX=',TRANSMAX
         IF(CHRIGIDROTT.AND.CHRMMT) WRITE(*,'(A,F8.4)')'  ROTMAX=',ROTMAX
         IF(.NOT.FIXTEMP(JP)) WRITE(*,'(A,I2,A,F8.4)') '[',JP,'] Temperature is now:',TEMP(JP)
         IF (RIGID) WRITE(*,'(A,I2,A,F12.6,A,F12.6)') 
     1           '[',JP,'] Maximum rigid body rotational move is now ',OSTEP(JP)
      ENDIF
      IF (FIXD) WRITE(*,'(A,I2,A,I4)') '[',JP,'] hard sphere collision moves=',NHSMOVE
C
      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 
C
      RETURN
      END

C
C
      SUBROUTINE REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
      USE commons
      IMPLICIT NONE
      LOGICAL STAY, REJECT, METROPOLIS
      INTEGER J1, J2, NREN, ITERATIONS, NQTOT, NLAST, JACCPREV, NSTEPREN, J3, BRUN, QDONE
      DOUBLE PRECISION POTEL, RMIN, RCOORDS(3*NATOMS), RVAT(NATOMS), RANDOM, DPRAND, RMINO, RCOORDSO(3*NATOMS), RVATO(NATOMS),
     1                 TIME, XIP, DUMMY, SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

      STAY=.FALSE.
      IF (POTEL.LT.RMIN) THEN
         RMIN=POTEL          
         DO J2=1,3*NATOMS
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            RVAT(J2)=VAT(J2,1)
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(*,'(A,2G20.10)') 'RMIN,POTEL=',RMIN,POTEL
C     PRINT*,'J1,JACCPREV+NRENSTUCK,NLAST+NREN=',J1,JACCPREV+NRENSTUCK,NLAST+NREN
      IF ((J1.GE.JACCPREV+NRENSTUCK).OR.(J1.GE.NLAST+NREN)) THEN
         JACCPREV=J1
         NLAST=J1
         RANDOM=DPRAND()
         METROPOLIS=DEXP(-(RMIN-RMINO)/TRENORM).GT.RANDOM
         REJECT=.FALSE.
C
C  Taboo list for renormalised energies. Skip if the step is going to be rejected by
C  the Metropolis condition.
C
         IF (TABOOT.AND.METROPOLIS) THEN
            IF (NSTEPREN.EQ.0) NT(1)=0
            CALL NEWINERTIA(RCOORDS,NATOMS,NATOMS,XIP)
            DO J1=1,NT(1)
               IF (DABS(RMIN-ESAVE(J1,1)).LT.ECONV) THEN
                  IF (DABS(XIP-XINSAVE(J1,1))/(XIP+XINSAVE(J1,1)).LT.1.0D-2) THEN
                     REJECT=.TRUE.
                     GOTO 20
                  ELSE
                     PRINT*,'Energies nearly degenerate:',RMIN,ESAVE(J1,1)
                     PRINT*,'But  different  structures:',XIP,XINSAVE(J1,1)
                  ENDIF
               ENDIF
               IF (RMIN.LT.ESAVE(J1,1)) THEN
                  NT(1)=MIN(NT(1)+1,NTAB)
                  DO J3=NT(1),J1+1,-1
                     ESAVE(J3,1)=ESAVE(J3-1,1)
                     XINSAVE(J3,1)=XINSAVE(J3-1,1)
                  ENDDO
                  ESAVE(J1,1)=RMIN
                  XINSAVE(J1,1)=XIP
                  GOTO 20
               ENDIF
            ENDDO
            
            NT(1)=MIN(NT(1)+1,NTAB)
            ESAVE(NT(1),1)=RMIN
            XINSAVE(NT(1),1)=XIP

20          CONTINUE

            WRITE(*,'(A,I10)') ' Number of entries in taboo list=',NT(1)
            IF (DEBUG) THEN
               WRITE(*,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
C
C  Accept/reject for renormalised step
C
         IF (METROPOLIS.AND.(.NOT.REJECT)) THEN
            IF (NSTEPREN.GT.0) WRITE(*,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' accepted'
            NREN=MAX(NREN/1.1D0,NRENORM/2.0D0)
            RMINO=RMIN
            IF (.NOT.STAY) THEN
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=RCOORDS(J2)
                  COORDS(J2,1)=RCOORDS(J2)
               ENDDO
               DO J2=1,NATOMS
                  RVATO(J2)=RVAT(J2)
               ENDDO
            ELSE
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=COORDSO(J2,1)
                  COORDS(J2,1)=COORDSO(J2,1)
               ENDDO
            ENDIF
         ELSE
            IF (REJECT) THEN
               WRITE(*,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected by taboo criterion'
            ELSE
               WRITE(*,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected'
            ENDIF
            DO J2=1,3*NATOMS
               COORDS(J2,1)=RCOORDSO(J2)
            ENDDO
            DO J2=1,NATOMS
               RVAT(J2)=RVATO(J2)
            ENDDO
            NREN=NREN*1.1D0
         ENDIF
         NSTEPREN=NSTEPREN+1
         IF (NSTEPREN.EQ.1) WRITE(*,'(A,G20.10)') ' first renorm energy is ',RMIN
         WRITE(*,'(A,I6)') ' renormalisation interval is now ',NREN
C
C  Longer renorm step
C
10       IF ((XMOVERENORM.GT.3.0D0).OR.FIXD) THEN
            CALL HSMOVE(COORDS,1,INT(XMOVERENORM))
         ELSE
            DUMMY=STEP(1)
            STEP(1)=XMOVERENORM
            CALL TAKESTEP(1)
            STEP(1)=DUMMY
         ENDIF
         CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
C
C  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
C
         IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=COORDSO(J2,1)
            ENDDO
            GOTO 10
         ENDIF
         NQTOT=NQTOT+1
         WRITE(*,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Renorm Qu ',NQ(1),' E=',
     1        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
         DO J2=1,3*NATOMS
            COORDSO(J2,1)=COORDS(J2,1)
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            VATO(J2,1)=VAT(J2,1)
            RVAT(J2)=VAT(J2,1)
         ENDDO
         EPREV(1)=POTEL
         RMIN=POTEL
      ENDIF

      RETURN
      END
C
C  Reseed if the energy has not improved by more than ECONV over the
C  last NRELAX mc steps.
C  If AVOID is true then save the energy and coordinates of the lowest
C  minimum achieved before each restart and restart if we get too close
C  to any one of them using bipartite matching and mind. Note that bipartite
C  matching and mind can give a local minimum of distance if the optimal 
C  permutation-inversion isn't found. Using ORIENT to effect a standard
C  orientation first seems to help. It should at least ensure that permutation-inversion
C  isomers are always found. 
C  Would it be possible to just use the energy and inertia components to identify
C  permutation-inversion isomers instead of MINPERM and MIND?
C  This would be like using a small value of AVOIDDIST, which doesn't seem to be
C  as good.
C
      SUBROUTINE NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J3, PERM(NATOMS), NPERM
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2,
     1                 RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS), DPRAND, FCOORDS(3*NATOMS), XMSBSAVE(3*NATOMS),
     2                 DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS), EBESTP, BESTCOORDSP(3*NATOMS), WORSTRAD
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), NORBIT1, NORBIT2, INVERT, NPOSITION
      LOGICAL RES1, RES2

      SR3=DSQRT(3.0D0)
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
C
C  Reseed if the energy has not improved in the last NRELAX mc cycles,
C  or if the current minimum is too close to one of the NMSBSAVE structures
C  saved in MSBCOORDS.
C
C  Instead of using the current minimum, employ the current minimum in 
C  BESTCOORDS for the AVOID check. Then we need only do the AVOID check when we have
C  a new minimum in BESTCOORDS, i.e. J1.EQ.JBEST(JP).
C
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     PRINT '(A,I5,2G17.7,3I5,L8)','J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
C     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
         FCOORDS(1:3*NATOMS)=DUMMY(1:3*NATOMS)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
C
C  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
C  but the coordinates in XMSB do. DISTANCE is the distance in this case.
C
! If the energy is lower no reseeding regardless of separation ? DJW
!           PRINT '(A,2G20.10,L8)','POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV=',
!    1                              POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            INVERT=1
30          XMSB(1:3*NATOMS)=INVERT*MSBCOORDS(1:3*NATOMS,J2)
C
C  If INVERT=-1 we need to realign XMSB
C
            IF (INVERT.EQ.-1) THEN
               CALL MYORIENT(XMSB,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG) ! saved MSBCOORDS were already oriented
               XMSB(1:3*NATOMS)=DUMMY(1:3*NATOMS)
            ENDIF
            XMSBSAVE(1:3*NATOMS)=XMSB(1:3*NATOMS)
10          CALL MINPERM(NATOMS, FCOORDS, XMSB, BOXLX, BOXLY, BOXLZ, PERIODIC, PERM, DISTANCE, DIST2, WORSTRAD)
! 10          CALL BIPARTITE(NATOMS, FCOORDS, XMSB, PERM, DISTANCE, DIST2, WORSTRAD)
            DISTANCE=SQRT(DISTANCE)
            DUMMY(1:3*NATOMS)=XMSB(1:3*NATOMS)
            NPERM=0
            DO J3=1,NATOMS
               XMSB(3*(J3-1)+1)=DUMMY(3*(PERM(J3)-1)+1)
               XMSB(3*(J3-1)+2)=DUMMY(3*(PERM(J3)-1)+2)
               XMSB(3*(J3-1)+3)=DUMMY(3*(PERM(J3)-1)+3)
               IF (J3.NE.PERM(J3)) THEN
!                 IF (DEBUG) WRITE(*,'(A,I5,A,I5)') 'permute atoms ',J3,' and ',PERM(J3)
                  NPERM=NPERM+1
               ENDIF
            ENDDO
            IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,G20.10)') 'distance from structure ',INVERT*J2,
     1         ' after permuting ',NPERM,' pairs of atoms=',DISTANCE
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(*,'(A,I6,A,G20.10,A,F10.3)') 'Current minimum is too close to saved structure ',
     1             INVERT*J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
C
C  Distance minimisations with respect to Euler angles and centre-of-mass.
C  Coordinates in XMSB are reset by mind (second argument).
C  
            IF (NPERM.NE.0) THEN 
               CALL MINDGMIN(FCOORDS,XMSB,NATOMS,DISTANCE,PERIODIC,TWOD)
               IF (DEBUG) WRITE(*,'(A,G20.10)') 
     1         'distance after mind=                                                 ',DISTANCE
               IF (DISTANCE.LT.AVOIDDIST) THEN
                  RES2=.TRUE.
                  WRITE(*,'(A,I6,A,G20.10,A,F10.3)') 'Current minimum is too close to saved structure ',
     1                INVERT*J2,' with energy ',MSBE(J2),' dist=',DISTANCE
                  GOTO 20
               ENDIF
               GOTO 10
            ELSE
               IF (INVERT.EQ.1) THEN
C
C  Now try the enantiomer if necessary.
C
                  INVERT=-1
                  GOTO 30
               ENDIF
               IF (ABS(POTEL-MSBE(J2)).LT.ECONV) THEN
                  WRITE(*,'(A,I6,A,2G20.10)') 'WARNING current energy and energy of structure ',J2,
     1                ' are ',POTEL,MSBE(J2)
C                 PRINT*,NATOMS
C                 PRINT*,'current coords:'
C                 WRITE(*,'(3G20.10)') (COORDS(J3,JP),J3=1,3*NATOMS)
C                 PRINT*,NATOMS
C                 PRINT*,'FCOORDS:'
C                 WRITE(*,'(3G20.10)') (FCOORDS(J3),J3=1,3*NATOMS)
C                 PRINT*,NATOMS
C                 PRINT*,'MSB coords:'
C                 WRITE(*,'(3G20.10)') (MSBCOORDS(J3,J2),J3=1,3*NATOMS)
C                 PRINT*,NATOMS
C                 PRINT*,'XMSBSAVE coords:'
C                 WRITE(*,'(3G20.10)') (XMSBSAVE(J3),J3=1,3*NATOMS)
C                 STOP
               ENDIF
            ENDIF
         ENDDO savedloop
      ENDIF
20    CONTINUE
      IF (RES1.OR.RES2) THEN
C
C  Does not seem necessary to save MSB data for restarts.
C  If we are reseeding because RES2 is true then:
C  (1) If the AVOID test is based on the distance for BESTCOORDS, add
C      these to the  MSB list if the distance is > 0.01D0, i.e. a different structure.
C      Otherwise, add the coordinates of the previous lowest energy minimum to
C      the MSB list; these have been saved in BESTCOORDSP and the energy in EBESTP.
C  (2) Alternatively, if the AVOID test is based on the coordinates of the current
C      structure, we should add these to the MSB list instead. Testing
C      every step involves an overhead, which grows with the number of saved taboo
C      structures. We could make the taboo list cyclic, so that new structures to be
C      avoided replace old ones at the beginning of the list when the length of the
C      list is exceeded. We are currently only doing the bipartite matching test if
C      the energy of the lowest minimum since the last reseeding has changed, which
C      saves time.
C
C  Change to a cyclic AVOID list 30/12/04.
C
         IF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.01D0))) THEN ! new condition
!!!!!!      IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
!!             IF (NMSBSAVE.GT.MAXSAVE) NMSBSAVE=1 ! cycle
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
C              PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
               MSBE(NPOSITION)=EBEST(JP)
               FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
!              MSBCOORDS(1:3*NATOMS,NPOSITION)=BESTCOORDS(1:3*NATOMS,JP)
               WRITE(*,'(A,I6,A,G20.10)') 'Moving current best minimum to position ',NPOSITION,
     1                             ' in the AVOID list E=',EBEST(JP)
!!!!!!      ENDIF
C           OPEN(UNIT=34,FILE='MSBdata',POSITION='APPEND')
C           WRITE(34,'(G20.10)') MSBE(NPOSITION) 
C           WRITE(34,'(3G20.10)') BESTCOORDS(1:3*NATOMS,JP)
C           CLOSE(34)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
!!!!!!      IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
!!             IF (NMSBSAVE.GT.MAXSAVE) NMSBSAVE=1 ! cycle
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
C              PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
               MSBE(NPOSITION)=EBESTP
C              FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! BUG ??
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
!              MSBCOORDS(1:3*NATOMS,NPOSITION)=BESTCOORDSP(1:3*NATOMS)
               WRITE(*,'(A,I6,A,G20.10)') 'Moving previous best minimum to position ',NPOSITION,
     1                             ' in the AVOID list E=',EBESTP
!!!!!!      ENDIF
         ENDIF ! end new condition
         IF (NHSRESTART.GT.0) THEN
            IF (RES1) WRITE(*,'(A,I8,A)') 'Energy has not improved since step ',JBEST(JP),' - perturbing'
            IF (RES2) WRITE(*,'(A,I8,A)') 'Reseeding due to taboo condition'
            CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
         ELSE
            IF (RES1) WRITE(*,'(A,I8,A)') 'Energy has not improved since step ',JBEST(JP),' - reseeding'
            IF (RES2) WRITE(*,'(A,I8,A)') 'Reseeding due to taboo condition'
            DO J2=1,3*NATOMS
               RANDOM=(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
            ENDDO
         ENDIF
C
C  This step will be accepted because EPREV(JP)=POTEL, so we should
C  not need to save COORDSO and VATO. Should reset EBEST and BESTCOORDS, though.
C
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=J1
         EPREV(JP)=POTEL ! +1.0D0
         EPPREV(JP)=0.0D0
         NSYMREM=0
      ENDIF

      RETURN
      END

C
C  Reseed if the energy has not improved by more than ECONV over the
C  last NRELAX mc steps.
C  If AVOID is true then save the energy and coordinates of the lowest
C  minimum achieved before each restart and restart if we get too close
C  to any one of them using bipartite matching and mind. Note that bipartite
C  matching and mind can give a local minimum of distance if the optimal 
C  permutation-inversion isn't found. Using ORIENT to effect a standard
C  orientation first seems to help. It should at least ensure that permutation-inversion
C  isomers are always found. 
C
C  In NEWRES2 do the AVOID check for every call.
C
      SUBROUTINE NEWRES2(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J3, PERM(NATOMS), NPERM
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2,
     1                 RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS), DPRAND, FCOORDS(3*NATOMS),
     2                 DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS), EBESTP, BESTCOORDSP(3*NATOMS), WORSTRAD
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), NORBIT1, NORBIT2, INVERT, NPOSITION
      LOGICAL RES1, RES2

      SR3=DSQRT(3.0D0)
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
C
C  Reseed if the energy has not improved in the last NRELAX mc cycles,
C  or if the current minimum is too close to one of the NMSBSAVE structures
C  saved in MSBCOORDS.
C
C  Employ the current minimum in the AVOID check.
C
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     PRINT '(A,I5,2G17.7,3I5,L8)','J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
C     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
!     IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
      IF ((.NOT.RES1).AND.AVOID.AND.(NMSBSAVE.GT.0)) THEN 
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
         FCOORDS(1:3*NATOMS)=DUMMY(1:3*NATOMS)
!        PRINT*,'in endif NMSBSAVE=',NMSBSAVE
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
C
C  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
C  but the coordinates in XMSB do. DISTANCE is the distance in this case.
C
! If the energy is lower no reseeding regardless of separation ? DJW
C
!           PRINT '(A,I5,2G20.10)','in savedloop J2,POTEL,MSBE=',J2,POTEL,MSBE(J2)
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            INVERT=1
30          XMSB(1:3*NATOMS)=INVERT*MSBCOORDS(1:3*NATOMS,J2)
C
C  If INVERT=-1 we need to realign XMSB
C
            IF (INVERT.EQ.-1) THEN
               CALL MYORIENT(XMSB,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG) ! saved MSBCOORDS were already oriented
               XMSB(1:3*NATOMS)=DUMMY(1:3*NATOMS)
            ENDIF
10          CALL MINPERM(NATOMS, FCOORDS, XMSB, BOXLX, BOXLY, BOXLZ, PERIODIC, PERM, DISTANCE, DIST2, WORSTRAD)
! 10          CALL BIPARTITE(NATOMS, FCOORDS, XMSB, PERM, DISTANCE, DIST2, WORSTRAD)
            DISTANCE=SQRT(DISTANCE)
!           PRINT '(A,I5,3G20.10)','J2,MSBE,DISTANCE,DIST2=',J2,MSBE(J2),DISTANCE,DIST2
            DUMMY(1:3*NATOMS)=XMSB(1:3*NATOMS)
            NPERM=0
            DO J3=1,NATOMS
               XMSB(3*(J3-1)+1)=DUMMY(3*(PERM(J3)-1)+1)
               XMSB(3*(J3-1)+2)=DUMMY(3*(PERM(J3)-1)+2)
               XMSB(3*(J3-1)+3)=DUMMY(3*(PERM(J3)-1)+3)
               IF (J3.NE.PERM(J3)) THEN
!                 IF (DEBUG) WRITE(*,'(A,I5,A,I5)') 'permute atoms ',J3,' and ',PERM(J3)
                  NPERM=NPERM+1
               ENDIF
            ENDDO
            IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,G20.10)') 'distance from structure ',INVERT*J2,
     1         ' after permuting ',NPERM,' pairs of atoms=',DISTANCE
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(*,'(A,I6,A,G20.10,A,F10.3)') 'Current minimum is too close to saved structure ',
     1             INVERT*J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
C
C  Distance minimisations with respect to Euler angles and centre-of-mass.
C  Coordinates in XMSB are reset by mind (second argument).
C  
            IF (NPERM.NE.0) THEN 
               CALL MINDGMIN(FCOORDS,XMSB,NATOMS,DISTANCE,PERIODIC,TWOD)
               IF (DEBUG) WRITE(*,'(A,G20.10)') 
     1         'distance after mind=                                                 ',DISTANCE
               IF (DISTANCE.LT.AVOIDDIST) THEN
                  RES2=.TRUE.
                  WRITE(*,'(A,I6,A,G20.10,A,F10.3)') 'Current minimum is too close to saved structure ',
     1                INVERT*J2,' with energy ',MSBE(J2),' dist=',DISTANCE
                  GOTO 20
               ENDIF
               GOTO 10
            ELSE
               IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,G20.10)') 'distance from structure ',INVERT*J2,
     1         ' after permuting ',NPERM,' pairs of atoms=',DISTANCE
               IF (INVERT.EQ.1) THEN
C
C  Now try the enantiomer if necessary.
C
                  INVERT=-1
                  GOTO 30
               ENDIF
               IF (ABS(POTEL-MSBE(J2)).LT.ECONV) THEN
                  WRITE(*,'(A,I6,A,2G20.10)') 'WARNING current energy and energy of structure ',J2,
     1                ' are ',POTEL,MSBE(J2)
C                 PRINT*,'current coords:'
C                 WRITE(*,'(3G20.10)') (COORDS(J3,JP),J3=1,3*NATOMS)
C                 PRINT*,'FCOORDS:'
C                 WRITE(*,'(3G20.10)') (FCOORDS(J3),J3=1,3*NATOMS)
C                 PRINT*,'MSB coords:'
C                 WRITE(*,'(3G20.10)') (MSBCOORDS(J3,J2),J3=1,3*NATOMS)
C                 PRINT*,'XMSB coords:'
C                 WRITE(*,'(3G20.10)') (XMSB(J3),J3=1,3*NATOMS)
C                 STOP
               ENDIF
            ENDIF
         ENDDO savedloop
      ENDIF
20    CONTINUE
      IF (RES1.OR.RES2) THEN
C
C  If we are reseeding because RES2 is true then add the coordinates of the current
C  minimum and the current MSB to the AVOID list.
C
C  Change to a cyclic AVOID list 30/12/04.
C
         IF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.01D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
!           PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)
            CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
            WRITE(*,'(A,I6,A,G20.10)') 'Moving current best minimum to position ',NPOSITION,
     1                          ' in the AVOID list E=',EBEST(JP)
            IF (RES2.AND.(J1.NE.JBEST(JP))) THEN ! add current minimum as well as MSB
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
!              PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
               MSBE(NPOSITION)=POTEL
               FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(*,'(A,I6,A,G20.10)') 'Moving current minimum to position ',NPOSITION,
     1                          ' in the AVOID list E=     ',POTEL
            ENDIF
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
            IF (J1.EQ.JBEST(JP)) THEN ! current minimum is current best
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
!              PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
               MSBE(NPOSITION)=EBESTP
!              FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ??? BUG
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(*,'(A,I6,A,G20.10)') 'Moving previous best minimum to position ',NPOSITION,
     1                                ' in the AVOID list E=',EBESTP
            ELSE
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
!              PRINT*,'NMSBSAVE,MAXSAVE,NPOSITION=',NMSBSAVE,MAXSAVE,NPOSITION
               MSBE(NPOSITION)=EBEST(JP)
               FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(*,'(A,I6,A,G20.10)') 'Moving current best minimum to position ',NPOSITION,
     1                          ' in the AVOID list E=',EBEST(JP)
            ENDIF
         ENDIF ! end new condition
         IF (NHSRESTART.GT.0) THEN
            IF (RES1) WRITE(*,'(A,I8,A)') 'Energy has not improved since step ',JBEST(JP),' - perturbing'
            IF (RES2) WRITE(*,'(A,I8,A)') 'Reseeding due to taboo condition'
            CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
         ELSE
            IF (RES1) WRITE(*,'(A,I8,A)') 'Energy has not improved since step ',JBEST(JP),' - reseeding'
            IF (RES2) WRITE(*,'(A,I8,A)') 'Reseeding due to taboo condition'
            DO J2=1,3*NATOMS
               RANDOM=(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
            ENDDO
         ENDIF
C
C  This step will be accepted because EPREV(JP)=POTEL, so we should
C  not need to save COORDSO and VATO. Should reset EBEST and BESTCOORDS, though.
C
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=J1
         EPREV(JP)=POTEL ! +1.0D0
         EPPREV(JP)=0.0D0
         NSYMREM=0
      ENDIF

      RETURN
      END

