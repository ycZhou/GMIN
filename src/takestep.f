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
      SUBROUTINE TAKESTEP(NP)
      USE commons
      IMPLICIT NONE
      

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS),
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS),
     2                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, NP, J3, JMAX2, RANATOM, INDEXD(NATOMS)
C     LOGICAL TMOVE(NPAR), OMOVE(NPAR)
C     COMMON /MOVE/ TMOVE, OMOVE
C     INTEGER NQTOT
C     COMMON /TOT/ NQTOT

C
C  This call can be used to keep the random numbers the same for parallel
C  runs - only for testing!
C
C     CALL SDPRND(NQTOT)

C     IF (NSYMREM.GT.0) THEN
C        OPEN(UNIT=77,FILE='coords.latest.xyz',STATUS='UNKNOWN')
C        WRITE(77,'(I6)') NATOMS
C        WRITE(77,'(A)') ' '
C        WRITE(77,'(A2,2X,3G20.10)') ('LA',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),COORDS(3*(J2-1)+3,NP),J2=1,NATOMS)
C        CLOSE(77)
C     ENDIF
C
C  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
C  outside the permitted radius. Then it may be impossible to take a step that keeps all the
C  atoms inside.
C
C     PRINT*,'in takestep VAT'
C     PRINT '(6G15.5)',VAT(1:NATOMS,NP)
C     PRINT*,'in takestep COORDS'
C     PRINT '(6G15.5)',COORDS(1:NATOMS,NP)

      IF ((.NOT.NORESET).AND.(.NOT.PERMUTE).AND.(.NOT.DIFFRACTT).AND.(.NOT.BLNT).AND.(.NOT.PERIODIC)
     &     .AND.(.NOT.GAUSST)) THEN
         DO J1=1,NATOMS
            IF ((.NOT.RIGID).OR.(J1.LE.NATOMS/2)) THEN
               J2=3*J1
               DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
               IF (DUMMY2.GT.RADIUS) THEN
                  WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                  PRINT*,'initial coordinate outside container - increase container radius'
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      DO J1=1,3*(NATOMS-NSEED)
         COORDSO(J1,NP)=COORDS(J1,NP)
      ENDDO
      DO J1=1,NATOMS
         VATO(J1,NP)=VAT(J1,NP)
      ENDDO
      IF (WENZEL) THEN
11       RANDOM=(DPRAND()-0.5D0)*2.0D0
         COORDS(1,NP)=COORDSO(1,NP)+STEP(NP)*RANDOM
         IF ((COORDS(1,NP).GT.1.0D0).OR.(COORDS(1,NP).LT.0.0D0)) GOTO 11
12       RANDOM=(DPRAND()-0.5D0)*2.0D0
         COORDS(2,NP)=COORDSO(2,NP)+STEP(NP)*RANDOM
         IF ((COORDS(2,NP).GT.1.0D0).OR.(COORDS(2,NP).LT.0.0D0)) GOTO 12
         RETURN
      ELSE IF (PERMUTE) THEN
         RETURN
      ENDIF
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS
      YMASS=YMASS/NATOMS
      ZMASS=ZMASS/NATOMS
C
C  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
C  and the pair energy of the most tightly bound atom, VMIN. An angular step is
C  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
C  DMAX (or CMMAX from CM of the cluster).
C
      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,NP).GT.VMAX) THEN
            VMAX=VAT(J1,NP)
            JMAX=J1
         ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
            VMAX2=VAT(J1,NP)
            JMAX2=J1
         ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO
C
C  If DECAY is true then select an atom at random, move this one randomly
C  by the maximum allowed amount, and move the others in random directions
C  of decaying magnitude, depending on how far they are from the chosen atom.
C
      IF (DECAY) THEN
9        RANATOM=NINT(0.5D0+NATOMS*DPRAND())
         IF (RANATOM.EQ.JMAX) GOTO 9 ! don't choose the atom that might undergo a surface move
         WRITE(*,'(A,I6)') 'atom undergoing maximum displacement is number ',RANATOM
         DO J1=1,NATOMS-NSEED
            DUMMY=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+
     1             (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2             (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
            EXPDECAY(J1)=EXP(-DECAYPARAM*DUMMY)
         ENDDO
      ENDIF
C
C  If COOP is true then select an atom at random, move this one randomly
C  by the maximum allowed amount, and move its NCOOP nearest neighbours by
C  the same displacement.
C
      IF (COOP) THEN
8        RANATOM=NINT(0.5D0+NATOMS*DPRAND())
         IF (RANATOM.EQ.JMAX) GOTO 8 ! don't choose the atom that might undergo a surface move
         IF (DEBUG) WRITE(*,'(A,I6)') 'randomly selected atom is number ',RANATOM
         DO J1=1,NATOMS-NSEED
            CDIST(J1)=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+
     1                 (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2                 (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
            INDEXD(J1)=J1
         ENDDO
         CALL SORT4(NATOMS-NSEED,NATOMS,CDIST,INDEXD)
         RANDOMX=(DPRAND()-0.5D0)*2.0D0
         RANDOMY=(DPRAND()-0.5D0)*2.0D0
         RANDOMZ=(DPRAND()-0.5D0)*2.0D0
         DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
         RANDOMX=RANDOMX/DUMMY2
         RANDOMY=RANDOMY/DUMMY2
         RANDOMZ=RANDOMZ/DUMMY2
      ENDIF

      DO J1=1,NATOMS-NSEED
C        IF (NMOVE.GT.0) THEN
C           IF (1.0D0*NMOVE/1.0D0*NATOMS.LT.DPRAND()) GOTO 13
C        ENDIF
10       J2=3*J1
         LOCALSTEP=STEP(NP)

         IF (RIGID.AND.(J1.GT.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
         ELSE IF (RIGID.AND.(J1.LE.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
         ENDIF
C
C  Angular move block.
C  If NORESET is .TRUE. then VAT won;t be set, so we should skip this block.
C
         IF ((((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX))).AND.(.NOT.RIGID).AND.(.NOT.BLNT).AND.
     &         (.NOT.DIFFRACTT).AND.(.NOT.GAUSST)
     &        .AND.(.NOT.NORESET).AND.(.NOT.PERIODIC).AND.(.NOT.THOMSONT)) THEN

           IF (DEBUG) WRITE(*,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 'angular move for atom ',J1,
     &           ' V=',VMAX,' Vmin=',VMIN,' next most weakly bound atom is ',JMAX2,' V=',VMAX2
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0

C          COORDS(J2-2,NP)=DMAX*DSIN(THETA)*DCOS(PHI)
C          COORDS(J2-1,NP)=DMAX*DSIN(THETA)*DSIN(PHI)
C          COORDS(J2,NP)=  DMAX*DCOS(THETA)
C          IF (DPRAND().GT.0.5D0) THEN

C
C  Evaporation is judged from the origin, not the centre of mass. We don't want the
C  angular move to cause evaporation. Obviously this will cause problems if we have a cluster that drifts
C  away from the origin up to the container radius.  
C
              COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
              COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
              COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
              IF (DUMMY.GT.RADIUS) THEN
                 DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
                 COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
                 COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
                 COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
              ENDIF
C          ELSE
C             COORDS(J2-2,NP)=XMASS+(COORDS(3*(JMAX2-1)+1,NP)-XMASS)*(CMMAX+1.0D0)/CMDIST(JMAX2)
C             COORDS(J2-1,NP)=YMASS+(COORDS(3*(JMAX2-1)+2,NP)-YMASS)*(CMMAX+1.0D0)/CMDIST(JMAX2)
C             COORDS(J2,NP)=  ZMASS+(COORDS(3*(JMAX2-1)+3,NP)-ZMASS)*(CMMAX+1.0D0)/CMDIST(JMAX2)
C          ENDIF

C          ELSE IF ((.NOT.ROUND).AND.(((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX2)).AND.(DZTEST))) THEN
C            IF (DEBUG) WRITE(*,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 
C     1                'angular move for atom ',J1,' V=',VMAX,' Vmin=',VMIN,' this is the next most weakly bound atom'
C            THETA=DPRAND()*PI
C            PHI=DPRAND()*PI*2.0D0
C
CC          COORDS(J2-2,NP)=DMAX*DSIN(THETA)*DCOS(PHI)
CC          COORDS(J2-1,NP)=DMAX*DSIN(THETA)*DSIN(PHI)
CC          COORDS(J2,NP)=  DMAX*DCOS(THETA)
C 
C            COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
C            COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
C            COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
C
C  End of angular move block.
C
         ELSE IF ((NATOMS-NSEED.EQ.1).AND.(NATOMS.GT.1)) THEN
           IF (DEBUG) WRITE(*,'(A,I4,A,F12.4,A,2F12.4)') 
     1                'angular move for atom ',J1,' V=',VAT(J1,NP),' Vmin, Vmax=',VMIN,VMAX
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
           COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
           COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
           DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           IF (DUMMY.GT.RADIUS) THEN
              DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
              COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
              COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
              COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
           ENDIF
         ELSE IF (DECAY) THEN
            RANDOMX=(DPRAND()-0.5D0)*2.0D0
            RANDOMY=(DPRAND()-0.5D0)*2.0D0
            RANDOMZ=(DPRAND()-0.5D0)*2.0D0
            DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
            COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX*EXPDECAY(J1)/DUMMY2
            COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY*EXPDECAY(J1)/DUMMY2
            COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ*EXPDECAY(J1)/DUMMY2
         ELSE IF (COOP) THEN
            IF ((J1.LE.NCOOP+1).AND.(CDIST(J1).LT.COOPCUT)) THEN
               J2=3*INDEXD(J1)
               COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX
               COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY
               COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ
            ELSE
               COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
            ENDIF
         ELSE IF (.NOT.FIXD) THEN
C
C  Tried changing to scale steps according to distance from CM. Maximum
C  allowed shift is linear with this distance. Worse.
C  Now try moving atoms according to how strongly bound they are.
C
           RANDOM=DPRAND()
C          WRITE(*,'(A,3F20.10)') 'EFAC,RANDOM,EXP=',EFAC,RANDOM,EXP(-EFAC*(VAT(J1,NP)-VMAX)/(VMIN-VMAX))
C          PRINT*,'VMIN,VMAX,EFAC=',VMIN,VMAX,EFAC
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-2) ! scale gauss steps by 1/GKSMALL
              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-1)
              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2)
              COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
C          ELSE IF (EXP(-EFAC*(VAT(J1,NP)-VMIN)/(VMAX-VMIN)).GT.RANDOM) THEN
           ELSE 
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=1.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1,NP)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
C             COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
C             COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
C             COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*DUMMY
           ENDIF
C
C Stop atoms leaving the container in this step
C
           IF ((.NOT.PERIODIC).AND.(.NOT.AMBER).AND.(.NOT.(RIGID.AND.((J1.GT.NATOMS/2)))).AND.(.NOT.BLNT)
     1                        .AND.(.NOT.DIFFRACTT).AND.(.NOT.THOMSONT).AND.(.NOT.GAUSST)) THEN
C          IF ((.NOT.PERIODIC).AND.(.NOT.AMBER).AND.(.NOT.(RIGID.AND.(LOCALSTEP.EQ.0.0D0))).AND.(.NOT.BLNT)) THEN
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
C
C  Simply rescaling the radius of an atom that leaves the container will bias the sampling
C  of configuration space.
C
              IF (DUMMY.GT.RADIUS) THEN
                 COORDS(J2-2,NP)=COORDSO(J2-2,NP)
                 COORDS(J2-1,NP)=COORDSO(J2-1,NP)
                 COORDS(J2,NP)=COORDSO(J2,NP)
                 DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
                 IF (DEBUG) WRITE(*,'(A,I5,3F20.10)') 'J1,DUMMY,RADIUS=',J1,DUMMY,RADIUS,DUMMY2
                 GOTO 10
              ENDIF
C                COORDS(J2-2,NP)=COORDS(J2-2,NP)*DSQRT(RADIUS/DUMMY)
C                COORDS(J2-1,NP)=COORDS(J2-1,NP)*DSQRT(RADIUS/DUMMY)
C                COORDS(J2,NP)=COORDS(J2,NP)*DSQRT(RADIUS/DUMMY)
C             ENDIF
           ENDIF
         ENDIF
         IF (TOSI.OR.WELCH) THEN
            DO J3=1,J1-1
               DUMMY=(COORDS(J2-2,NP)-COORDS(3*(J3-1)+1,NP))**2
     1              +(COORDS(J2-1,NP)-COORDS(3*(J3-1)+2,NP))**2
     2              +(COORDS(J2  ,NP)-COORDS(3*(J3-1)+3,NP))**2
               IF (DUMMY.LT.1.0D0) GOTO 10 
            ENDDO
         ENDIF
C13       CONTINUE
      ENDDO
      IF (FIXD) CALL HSMOVE(COORDS,NP,NHSMOVE)
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE(COORDS,NP)
C     PRINT*,'NSYMREM in takestep=',NSYMREM
      IF (NSYMREM.GT.0) THEN
         CALL KEEPSYM(NP)
C        OPEN(UNIT=77,FILE='coords.latest.xyz',STATUS='OLD',POSITION='APPEND')
C        WRITE(77,'(I6)') NATOMS
C        WRITE(77,'(A)') ' '
C        WRITE(77,'(A2,2X,3G20.10)') ('LA',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),COORDS(3*(J2-1)+3,NP),J2=1,NATOMS)
C        CLOSE(77)
      ENDIF

      RETURN
      END

      SUBROUTINE KEEPSYM(NP)
      USE commons
      IMPLICIT NONE
      INTEGER NP, J2, SYMINDEX(NATOMS), J3, J4
      DOUBLE PRECISION LCOORDS(3*NATOMS), NEWQ(3*NATOMS), SYMDELTA(3*NATOMS), DELTA(3*NATOMS),  SYMOP1(3,3)
      DOUBLE PRECISION DUMMY, STEPLENGTH, NEWSTEPLENGTH, ODIST1, ODIST2, DMIN
      LOGICAL ASSIGNED(NATOMS), BAD

      LCOORDS(1:3*NATOMS)=COORDSO(1:3*NATOMS,NP)
      DELTA(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)-COORDSO(1:3*NATOMS,NP)
C     STEPLENGTH=SUM(DELTA(1:3*NATOMS)**2)
      SYMDELTA(1:3*NATOMS)=DELTA(1:3*NATOMS)
C
C  New algorithm - choose the closest unclaimed atom in each case, so that
C  no tolerances are involved. 
C
      DO J2=1,NSYMREM
         BAD=.FALSE.
         SYMOP1(1:3,1:3)=SYMREM(J2,1:3,1:3)
         CALL MATMULV(NEWQ,LCOORDS,SYMOP1,NATOMS,3,3)
C          DO J3=1,NATOMS
C             SYMINDEX(J3)=0
C             matchloop: DO J4=1,NATOMS
C                DUMMY=SQRT((LCOORDS(3*(J4-1)+1)-NEWQ(3*(J3-1)+1))**2+ 
C      &                    (LCOORDS(3*(J4-1)+2)-NEWQ(3*(J3-1)+2))**2+ 
C      &                    (LCOORDS(3*(J4-1)+3)-NEWQ(3*(J3-1)+3))**2)
C C              PRINT '(A,2I5,G20.10)','J4,J4,DUMMY=',J3,J4,DUMMY
C                IF (DUMMY.LT.SYMTOL5) THEN
C                   SYMINDEX(J3)=J4
C                   CYCLE matchloop
C                ENDIF
C             ENDDO matchloop
C             IF (SYMINDEX(J3).EQ.0) THEN
C                PRINT '(A,I5,A)','warning, no matching atom for number ',J3,' in KEEPSYM - exiting'
C C              STOP
C                RETURN ! give up - this probably shouldn;t be a hard fail because there are cases where
C             ENDIF     ! it probably could occur legitimately
C          ENDDO
         
         ASSIGNED(1:NATOMS)=.FALSE.
         SYMINDEX(1:NATOMS)=0
         DO J3=1,NATOMS
            DMIN=1.0D100
            DO J4=1,NATOMS
               DUMMY=(LCOORDS(3*(J4-1)+1)-NEWQ(3*(J3-1)+1))**2+
     &               (LCOORDS(3*(J4-1)+2)-NEWQ(3*(J3-1)+2))**2+
     &               (LCOORDS(3*(J4-1)+3)-NEWQ(3*(J3-1)+3))**2
C              PRINT '(A,2I5,2G15.5)','J3,J4,DUMMY,DMIN=',J3,J4,DUMMY,DMIN
               IF (DUMMY.LT.DMIN) THEN
                  IF (ASSIGNED(J4)) THEN
C                    PRINT '(2(A,I5),A,F12.3)','WARNING closest atom ',J4,' to image atom ',J3,
C    &                                       ' already assigned, dist=',SQRT(DUMMY)
                  ELSE
                     IF (SYMINDEX(J3).GT.0) ASSIGNED(SYMINDEX(J3))=.FALSE.
                     SYMINDEX(J3)=J4
                     ASSIGNED(J4)=.TRUE.
                     DMIN=DUMMY
                  ENDIF
               ENDIF
            ENDDO
            IF (DEBUG.AND.(SQRT(DMIN).GT.SYMTOL5)) THEN
               PRINT '(A,I5,F12.3)','WARNING closest image to atom ',J3,' distance=',SQRT(DMIN)
               BAD=.TRUE.
            ENDIF
         ENDDO
!        IF (BAD) THEN
!           OPEN(UNIT=1,FILE='keepsym.xyz',STATUS='UNKNOWN')
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'LCOORDS'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',LCOORDS(3*(J4-1)+1),LCOORDS(3*(J4-1)+2),LCOORDS(3*(J4-1)+3)
!           ENDDO
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'NEWQ'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',NEWQ(3*(J4-1)+1),NEWQ(3*(J4-1)+2),NEWQ(3*(J4-1)+3)
!           ENDDO
!           CLOSE(1)
!        ENDIF

         CALL MATMULV(NEWQ,DELTA,SYMOP1,NATOMS,3,3)
         DO J3=1,NATOMS
            SYMDELTA(3*(J3-1)+1)=SYMDELTA(3*(J3-1)+1)+NEWQ(3*(SYMINDEX(J3)-1)+1)
            SYMDELTA(3*(J3-1)+2)=SYMDELTA(3*(J3-1)+2)+NEWQ(3*(SYMINDEX(J3)-1)+2)
            SYMDELTA(3*(J3-1)+3)=SYMDELTA(3*(J3-1)+3)+NEWQ(3*(SYMINDEX(J3)-1)+3)
         ENDDO
      ENDDO

C     NEWSTEPLENGTH=SUM(SYMDELTA(1:3*NATOMS)**2)
C     SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)*SQRT(STEPLENGTH/NEWSTEPLENGTH)
C
C  Maintain CofM distance in symmetry adapted step.
C
      SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)/(1+NSYMREM)
      DO J2=1,NATOMS
         ODIST1=SUM(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)**2)
         ODIST2=SUM((COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)+SYMDELTA(3*(J2-1)+1:3*(J2-1)+3))**2)
         IF (ODIST2.NE.0.0D0) COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)=(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)
     &        +SYMDELTA(3*(J2-1)+1:3*(J2-1)+3)) * SQRT(ODIST1/ODIST2)
      ENDDO
    
      RETURN
      END
