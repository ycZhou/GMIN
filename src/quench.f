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
C
C  Conjugate gradient driver. 
C  CFLAG convergence test
C  CTEST checks for changes in chirality for AMBER runs
C
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
      USE commons
      use qmodule
      use porfuncs
      IMPLICIT NONE

      INTEGER I, J1, NSQSTEPS, NP, IFLAG, ITER, NOPT, J2
      DOUBLE PRECISION P(3*NATOMS),POTEL,TIME,EREAL,RBCOORDS(18),TMPCOORDS(3*NATOMS), DIST
      LOGICAL QTEST, CFLAG, RES, COMPON, evapreject
      DOUBLE PRECISION  GRAD(3*NATOMS), DUMMY, DUM(3*NATOMS), DISTMIN, SSAVE, DIST2
C     DOUBLE PRECISION  WORK(60*NATOMS)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0

      CHARACTER(LEN=80) DSTRING
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      LOGICAL GUIDECHANGET, GUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET
      common /ev/ evapreject
      DOUBLE PRECISION QSTART, QFINISH
      COMMON /Q4C/ QSTART, QFINISH
C
C  Data for the screen saver.
C
      INTEGER BRUN, QDONE
C
C  Turn on guiding potentials. These get turned off in potential.F when
C  the RMS force is small enough.
C
      SSAVE=STEP(NP)

11    IF (WELCH) TOSI=.TRUE.
      IF (PACHECO) AXTELL=.FALSE.
      IF (CPMD) SCT=.TRUE.
      IF (ZETT1.OR.ZETT2) THEN
         MORSET=.TRUE.
         RHO=6.0D0
      ENDIF
      IF (NATBT.AND.GUPTAT) GUIDET=.TRUE.
      IF (NATBT.AND.GUIDET) GUPTAT=.TRUE.
      NOPT=3*NATOMS
      IF (WENZEL) NOPT=2
C
C  QTEST is set for the final quenches with tighter convergence criteria.
C
      IF (QTEST) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      QDONE=0
      DO I=1,3*NATOMS
         P(I)=COORDS(I,NP)
      ENDDO
C     IF (TIP) THEN
C        WRITE(40,'(I6)') (NATOMS/2)*3
C        WRITE(40,70) NP,NQ(NP), EREAL, RMS
C        DO J2=1,NATOMS/2
C           CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
C    1           P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
C           WRITE(40,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
C           WRITE(40,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
C           WRITE(40,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
C        ENDDO
C     ENDIF


      IF (COMPRESST.AND.(.NOT.QTEST)) THEN
         COMPON=.TRUE.
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         POTEL=EREAL
         IF (.NOT.CFLAG) WRITE(*,'(A,I7,A)') ' WARNING - compressed quench ',NQ(NP),'  did not converge'
         WRITE(*,'(A,I7,A,F20.10,A,I5,A,F15.7,A,I4,A,F12.2)') 'Comp Q ',NQ(NP),' energy=',
     1              POTEL,' steps=',ITER,' RMS force=',RMS
      ENDIF
      COMPON=.FALSE.

10    IF (PERMUTE) THEN
         CALL POTENTIAL(P,GRAD,EREAL,.FALSE.,.FALSE.)
         CFLAG=.TRUE.
         ITER=1
         RMS=0.0D0
      ELSE IF (DL_POLY) THEN
C
C  Need to make DL_POLY input file for current coordinates.
C
         OPEN (UNIT=91,FILE='CONFIG',STATUS='OLD')
         OPEN (UNIT=92,FILE='config',STATUS='UNKNOWN')
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         DO J1=1,NATOMS
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(3G20.10)') P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3)
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
         ENDDO
         CLOSE(91)
         CLOSE(92)
         CALL SYSTEM('cp CONFIG CONFIG.old; cp config CONFIG')
         CALL SYSTEM('DLPOLY.X > output.DL_POLY ; tail -9 STATIS > energy')
         OPEN (UNIT=91,FILE='energy',STATUS='OLD')
         READ(91,*) EREAL
         WRITE(*,'(A,G20.10)') 'energy=',EREAL
         CLOSE(91)
         OPEN(UNIT=91,FILE='REVCON',STATUS='OLD')
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
         NATOMS=0
13       READ(91,'(A1)',END=14) DUMMY
         NATOMS=NATOMS+1
         READ(91,*) P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
C        WRITE(*,'(3G20.10)') P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         GOTO 13
14       CONTINUE
         CLOSE(91)
         CFLAG=.TRUE.
C
C  Read the coordinates of the minimised geometry into vector P.
C
C     ELSE IF (BFGS .AND.(.NOT.QTEST)) THEN
      ELSE IF (BFGS) THEN
C        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
         CALL DFPMIN(MAXIT,P,3*NATOMS,GMAX,ITER,EREAL,CFLAG)
      ELSEIF (TNT) THEN
C        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
         PRINT '(A)','subroutine tn does not compile with NAG/PG'
         STOP
C        CALL TN(IFLAG,3*NATOMS,P,EREAL,GRAD,WORK,60*NATOMS,GMAX,ITER,MAXIT,CFLAG,DEBUG)
      ELSEIF (CONJG) THEN
         CALL CGMIN(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSEIF (RKMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSEIF (BSMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSE
C        CALL CGMIN(5,P,CFLAG,ITER,EREAL,NP)
         IF (CHRMMT.AND.INTMINT) THEN
            CALL MYLBFGS(NINTS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ELSE IF (THOMSONT) THEN
            DO J1=1,NATOMS
               DIST=SQRT(COORDS(3*(J1-1)+1,NP)**2+COORDS(3*(J1-1)+2,NP)**2+COORDS(3*(J1-1)+3,NP)**2)
               COORDS(3*(J1-1)+1,NP)=COORDS(3*(J1-1)+1,NP)/DIST
               COORDS(3*(J1-1)+2,NP)=COORDS(3*(J1-1)+2,NP)/DIST
               COORDS(3*(J1-1)+3,NP)=COORDS(3*(J1-1)+3,NP)/DIST
               P(2*(J1-1)+1)=ACOS(COORDS(3*(J1-1)+3,NP))
               IF (ABS(COORDS(3*(J1-1)+3,NP)-COS(P(2*(J1-1)+1))).GT.1.0D-10) THEN
                  PRINT '(A)','inconsistent conversion for z'
                  STOP
               ENDIF
               IF (COORDS(3*(J1-1)+1,NP).EQ.0.0D0) THEN
                  IF (COORDS(3*(J1-1)+2,NP).GT.0.0D0) THEN
                     P(2*(J1-1)+2)=HALFPI
                  ELSE 
                     P(2*(J1-1)+2)=-HALFPI
                  ENDIF
               ELSE IF (COORDS(3*(J1-1)+2,NP).EQ.0.0D0) THEN
                  IF (COORDS(3*(J1-1)+1,NP).GT.0.0D0) THEN
                     P(2*(J1-1)+2)=0.0D0
                  ELSE 
                     P(2*(J1-1)+2)=2*HALFPI
                  ENDIF
               ELSE
                  P(2*(J1-1)+2)=ATAN(COORDS(3*(J1-1)+2,NP)/COORDS(3*(J1-1)+1,NP))
               ENDIF
               IF (ABS(COORDS(3*(J1-1)+1,NP)-SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                  P(2*(J1-1)+2)=P(2*(J1-1)+2)+2*HALFPI
                  IF (ABS(COORDS(3*(J1-1)+1,NP)-SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                     PRINT '(A)','inconsistent conversion for x'
                     STOP
                  ENDIF
               ENDIF
               IF (ABS(COORDS(3*(J1-1)+2,NP)-SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                  P(2*(J1-1)+2)=-P(2*(J1-1)+2)
                  IF (ABS(COORDS(3*(J1-1)+2,NP)-SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2))).GT.1.0D-5) THEN
                     PRINT '(A)','inconsistent conversion for y'
                     PRINT '(A,3G20.10)','x,y,z:      ',COORDS(3*(J1-1)+1,NP),COORDS(3*(J1-1)+2,NP),COORDS(3*(J1-1)+3,NP)
                     PRINT '(A,3G20.10)','theta,phi: ',P(2*(J1-1)+1),P(2*(J1-1)+2)
                     PRINT '(A,3G20.10)','x,y,z calc: ',SIN(P(2*(J1-1)+1))*COS(P(2*(J1-1)+2)),
     &                                                  SIN(P(2*(J1-1)+1))*SIN(P(2*(J1-1)+2)),
     &                                                  COS(P(2*(J1-1)+1))
                     STOP
                  ENDIF
               ENDIF
            ENDDO
            CALL MYLBFGS(2*NATOMS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
            TMPCOORDS(1:2*NATOMS)=P(1:2*NATOMS) 
            DO J1=1,NATOMS
               P(3*(J1-1)+1)=SIN(TMPCOORDS(2*(J1-1)+1))*COS(TMPCOORDS(2*(J1-1)+2))
               P(3*(J1-1)+2)=SIN(TMPCOORDS(2*(J1-1)+1))*SIN(TMPCOORDS(2*(J1-1)+2))
               P(3*(J1-1)+3)=COS(TMPCOORDS(2*(J1-1)+1))
            ENDDO
         ELSE
            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ENDIF
         if (evapreject) return
      ENDIF
      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(*,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(*,'(A,I2,A,I7,A)') '[',NP,'] WARNING - Quench ',NQ(NP),'  did not converge'
            ELSE
               WRITE(*,'(A,I7,A)') ' WARNING - Quench ',NQ(NP),'  did not converge'
            ENDIF
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.
      IF (TABOOT.AND.(.NOT.QTEST).AND.(.NOT.RENORM)) THEN
         CALL TABOO(EREAL,POTEL,P,NP,RES)
         IF (RES) GOTO 10
      ENDIF

C     PRINT*,'Taboo lists:'
C     DO J1=1,NPAR
C        PRINT*,'Parallel run ',J1
C        WRITE(*,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
C     ENDDO
C     PRINT*,'Inertia lists:'
C     DO J1=1,NPAR
C        PRINT*,'Parallel run ',J1
C        WRITE(*,'(6F15.7)') (XINSAVE(J2,J1),J2=1,NT(J1))
C     ENDDO

      IF (SAVEQ) CALL GSAVEIT(EREAL,P,NP)
      IF (QDONE.EQ.0) THEN
         PRINT*,'quench did not converge from starting coodinates:'
         WRITE(*,'(3G20.10)') (COORDS(J1,NP),J1=1,3*NATOMS)
      ENDIF
C
C  If EPSSPHERE is non-zero we are presumably doing a calculation of the 
C  energy density of local minima. We need to know the minimum distance
C  between the starting point and the quenched minima.
C
      IF ((EPSSPHERE.NE.0.0D0).OR.HIST) THEN
         DO J1=1,3*NATOMS
            DUM(J1)=COORDS(J1,NP)
         ENDDO
C
C  DUM is returned in the closest orientation to P; P should not change.
C  This is nearly the same mind as OPTIM! To execute a random walk we must take 
C  another step and minimise until the distance between the starting point
C  and the quench minimum is less than EPSSPHERE.
C
         CALL MINDGMIN(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD)
      ENDIF
C
C  Deal with EPSSPHERE sampling.
C
      IF (EPSSPHERE.NE.0.0D0) THEN
         IF ((DISTMIN.GT.EPSSPHERE).OR.(ABS(EREAL-EPREV(NP)).LE.ECONV)) THEN
            WRITE(*,'(A,F12.5,A,4F14.5)') 'step ',STEP(NP),' EREAL, EPREV, DISTMIN, EPSSPHERE=',
     1                                     EREAL, EPREV(NP), DISTMIN, EPSSPHERE
            DO J1=1,3*NATOMS
               COORDS(J1,NP)=COORDSO(J1,NP)
            ENDDO
            CALL TAKESTEP(NP)
            PRINT*,'reseeding step, maximum displacement reset to ',STEP(NP)
            GOTO 11
         ELSE
            WRITE(*,'(A,2F20.10)') 'valid step, DISTMIN, EPSSPHERE=',DISTMIN, EPSSPHERE
         ENDIF
      ENDIF
C
C  If we are provided with target minimum coordinates in file coords.target then
C  calculate the minimum distances. May be useful for algorithm development.
C  If we get close, we don;t want to escape without a hit!
C
!     IF (ALLOCATED(TCOORDS)) THEN
!        DO J1=1,NTARGETS
!           TMPCOORDS(1:3*NATOMS)=TCOORDS(J1,1:3*NATOMS)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2)
!           PRINT '(A,I5,A,F15.3,A,F15.3,A,F20.10)','for target structure ',J1,' dist=',SQRT(DUMMY),' dist2=',DIST2,' V=',POTEL
!        ENDDO
!        DO J1=1,MIN(NMSBSAVE,MAXSAVE)
!           TMPCOORDS(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J1)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2)
!           PRINT '(A,I5,A,F15.3,A,F15.3,A,F20.10)','for taboo  structure ',J1,' dist=',SQRT(DUMMY),' dist2=',DIST2,' V=',POTEL
!        ENDDO
!     ENDIF
C
C  NORESET true does not set the configuration point to the quench geometry
C  A relaxed frozen core does not get saved, but the lowest minima are saved
C  by GSAVEIT.
C
      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      IF (Q4T) CALL ORDERQ4(NATOMS,P,QFINISH)
C
C  Calling CENTRE here without an evaporation check can put particles
C  outside the container, and make a valid step in takestep impossible.
C
C     PRINT*,'Calling centre from quench'
C     IF ((.NOT.FIELDT).AND.(.NOT.SEEDT).AND.CENT) CALL CENTRE(COORDS,NP)

      IF (DUMPT) THEN
         IF (ARNO) THEN
            WRITE(40,'(I4)') NATOMS+2
            WRITE(40,70) NP,NQ(NP),EREAL,RMS
            WRITE(40,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(40,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            IF (NS.NE.0) WRITE(40,65) (P(I),I=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TIP) THEN
            WRITE(39,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(40,'(I6)') (NATOMS/2)*3
            WRITE(40,70) NP,NQ(NP), EREAL, RMS
            DO J2=1,NATOMS/2
               CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
     1              P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               WRITE(40,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(40,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(40,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ELSE
            WRITE(39,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(40,'(I4)') NATOMS
            WRITE(40,70) NP,NQ(NP), EREAL, RMS
70          FORMAT(1X,'[',I2,'] QUENCH NUMBER ',I6,' final energy=',F20.10,' RMS force=',E20.10)
            WRITE(40,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.NE.0) WRITE(40,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
80          FORMAT(A2,3F20.10)
         ENDIF
      ENDIF

      IF (SQUEEZET) THEN
         IF ((EREAL.GT.0.0D0).AND.(SQUEEZED.LT.1.0D0)) THEN
            SQUEEZED=2.0D0-SQUEEZED
            NSQSTEPS=NQ(NP)
         ELSE
            NSQSTEPS=100000
         ENDIF
         DO J1=1,3*NVEC
            VEC(J1)=VEC(J1)*SQUEEZED
         ENDDO
         IF (NQ(NP).GT.2*NSQSTEPS) SQUEEZET=.FALSE.
      ENDIF
    
      IF ((NQ(NP).GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(*,'(I6,A,G20.10)') NSSTOP,' quenches completed, setting coordinates to the lowest minimum, E=',QMIN(1)
         DO J1=1,3*NATOMS
            COORDS(J1,NP)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      END
