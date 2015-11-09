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
      SUBROUTINE GSEED
      USE commons
      IMPLICIT NONE
      
      DOUBLE PRECISION XMASS, YMASS, ZMASS, DIST, DUMMY, DMAX, 
     1                 DMIN, AX, BX
      INTEGER J1, I, J2


      OPEN(UNIT=10,FILE='seed',STATUS='OLD')
      NSEED=0
      DO J1=1,NATOMS
         READ(10,*,END=10) COORDS(3*(NATOMS-J1)+1,1),COORDS(3*(NATOMS-J1)+2,1),COORDS(3*(NATOMS-J1)+3,1)
         NSEED=NSEED+1
      ENDDO
10    CLOSE(10)
      WRITE(*,'(A,I4,A)') 'Read core from file seed containing ',NSEED,' atoms'
      IF (FREEZECORE) THEN
         WRITE(*,'(A,I4,A)') 'Core will be fixed during the first ',NSSTOP,' quenches'
      ELSE
         WRITE(*,'(A,I4,A)') 'Core will be relaxed but reset for the first ',NSSTOP,' quenches'
      ENDIF
C
C  Centre the seed.
C
      IF (CENT) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=NATOMS,NATOMS-NSEED+1,-1
            XMASS=XMASS+COORDS(3*(I-1)+1,1)
            YMASS=YMASS+COORDS(3*(I-1)+2,1)
            ZMASS=ZMASS+COORDS(3*(I-1)+3,1)
         ENDDO
         XMASS=XMASS/NSEED
         YMASS=YMASS/NSEED
         ZMASS=ZMASS/NSEED
         DO I=NATOMS,NATOMS-NSEED+1,-1
            COORDS(3*(I-1)+1,1)=COORDS(3*(I-1)+1,1)-XMASS
            COORDS(3*(I-1)+2,1)=COORDS(3*(I-1)+2,1)-YMASS
            COORDS(3*(I-1)+3,1)=COORDS(3*(I-1)+3,1)-ZMASS
         ENDDO
C
C  Centre the other atoms.
C
C        IF (NATOMS-NSEED.GT.1) THEN
C           XMASS=0.0D0
C           YMASS=0.0D0
C           ZMASS=0.0D0
C           DO I=1,NATOMS-NSEED
C              XMASS=XMASS+COORDS(3*(I-1)+1,1)
C              YMASS=YMASS+COORDS(3*(I-1)+2,1)
C              ZMASS=ZMASS+COORDS(3*(I-1)+3,1)
C           ENDDO
C           PRINT*,'NATOMS-NSEED=',NATOMS-NSEED
C           XMASS=XMASS/(NATOMS-NSEED)
C           YMASS=YMASS/(NATOMS-NSEED)
C           ZMASS=ZMASS/(NATOMS-NSEED)
C           DO I=1,NATOMS-NSEED
C              COORDS(3*(I-1)+1,1)=COORDS(3*(I-1)+1,1)-XMASS
C              COORDS(3*(I-1)+2,1)=COORDS(3*(I-1)+2,1)-YMASS
C              COORDS(3*(I-1)+3,1)=COORDS(3*(I-1)+3,1)-ZMASS
C           ENDDO
C        ENDIF
      ENDIF
C
C  Find the largest radius vector of the seed.
C
      DIST=0.0D0
      DO J1=NATOMS,NATOMS-NSEED+1,-1
         DUMMY=COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2
         IF (DUMMY.GT.DIST) DIST=DUMMY
      ENDDO
      DIST=1.2D0*DSQRT(DIST)
C
C  Shift the coordinates of the non-core atoms outside the core.
C
      DMAX=0.0D0
      DMIN=1.0D20
      DO J1=1,NATOMS-NSEED
         DUMMY=DSQRT(COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2)
         IF (DUMMY.GT.DMAX) DMAX=DUMMY
         IF (DUMMY.LT.DMIN) DMIN=DUMMY
      ENDDO
      IF (DABS(DMAX-DMIN).GT.1.0D-5) THEN
         AX=DIST+DMIN*(DIST-DMIN-DMAX)/(DMAX-DMIN)
         BX=(DMIN+DMAX-DIST)/(DMAX-DMIN)
      ELSE
         AX=DIST
         BX=0.0D0
      ENDIF
      DO J1=1,NATOMS-NSEED
         DUMMY=DSQRT(COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2)
         COORDS(3*J1-2,1)=COORDS(3*J1-2,1)*(AX+BX*DUMMY)/DUMMY
         COORDS(3*J1-1,1)=COORDS(3*J1-1,1)*(AX+BX*DUMMY)/DUMMY
         COORDS(3*J1,1)  =COORDS(3*J1,1)  *(AX+BX*DUMMY)/DUMMY
      ENDDO
      WRITE(*,75)
75    FORMAT('Coordinates:')
      WRITE(*,80) (COORDS(J1,1),J1=1,3*NATOMS)
80    FORMAT(3F15.5)
      IF (DUMPT) THEN
         WRITE(40,*) NATOMS
         WRITE(40,*) ' Initial coordinates'
         WRITE(40,45) (COORDS(J2,1),J2=1,3*(NATOMS-NSEED))
45       FORMAT('LA ',3F20.10)
         WRITE(40,46) (COORDS(J2,1),J2=3*(NATOMS-NSEED)+1,3*NATOMS)
46       FORMAT('LB',3F20.10)
      ENDIF

      DO J1=1,3*NATOMS
         COORDSO(J1,1)=COORDS(J1,1)
      ENDDO
      
      NS=NSEED

      RETURN
      END
