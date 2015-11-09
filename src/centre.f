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
C**************************************************************************
C
C  Subroutine CENTRE moves the centre of mass to the origin.
C
C*********************************************************************
C
      SUBROUTINE CENTRE(X,NP)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION XMASS, YMASS, ZMASS, DIST, DISTMAX, X(3*NATOMS,NPAR)
      INTEGER I, J1, NP

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      IF (RIGID) THEN
         DO I=1,NATOMS/2
            XMASS=XMASS+X(3*(I-1)+1,NP)
            YMASS=YMASS+X(3*(I-1)+2,NP)
            ZMASS=ZMASS+X(3*(I-1)+3,NP)
         ENDDO
         XMASS=2*XMASS/NATOMS
         YMASS=2*YMASS/NATOMS
         ZMASS=2*ZMASS/NATOMS
         DO I=1,NATOMS/2
            X(3*(I-1)+1,NP)=X(3*(I-1)+1,NP)-XMASS
            X(3*(I-1)+2,NP)=X(3*(I-1)+2,NP)-YMASS
            X(3*(I-1)+3,NP)=X(3*(I-1)+3,NP)-ZMASS
         ENDDO
      ELSE
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1,NP)
            YMASS=YMASS+X(3*(I-1)+2,NP)
            ZMASS=ZMASS+X(3*(I-1)+3,NP)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
C        PRINT*,'initial coordinates in centre:'
C        WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1,NP),X(3*(I-1)+2,NP),X(3*(I-1)+3,NP),I=1,NATOMS)
         DO I=1,NATOMS
            X(3*(I-1)+1,NP)=X(3*(I-1)+1,NP)-XMASS
            X(3*(I-1)+2,NP)=X(3*(I-1)+2,NP)-YMASS
            X(3*(I-1)+3,NP)=X(3*(I-1)+3,NP)-ZMASS
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(*,'(A,3F12.4)') 'centre of mass reset to the origin from ',XMASS,YMASS,ZMASS
C     PRINT*,'final coordinates in centre:'
C     WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1,NP),X(3*(I-1)+2,NP),X(3*(I-1)+3,NP),I=1,NATOMS)
C     IF (RIGID) RETURN
C
C  Check that all the atoms are in the container. If not then rescale.
C  Must not do this - it could change the system from a minimum!
C
C     DISTMAX=0.0D0
C     DO J1=1,NATOMS
C        DIST=X(3*J1-2,NP)**2+X(3*J1-1,NP)**2+X(3*J1,NP)**2
C        IF (DIST.GT.DISTMAX) DISTMAX=DIST
C     ENDDO
C     IF (DISTMAX.GT.RADIUS) THEN
C        PRINT*,'DISTMAX,RADIUS=',DISTMAX,RADIUS
C        DISTMAX=DSQRT(DISTMAX/RADIUS)*0.99D0
C        DO J1=1,NATOMS
C           PRINT*,'J1,dist=',J1,X(3*J1-2,NP)**2+X(3*J1-1,NP)**2+X(3*J1,NP)**2
C           X(3*J1-2,NP)=X(3*J1-2,NP)*DISTMAX
C           X(3*J1-1,NP)=X(3*J1-1,NP)*DISTMAX
C           X(3*J1,NP)  =X(3*J1,NP)  *DISTMAX
C        ENDDO
C     ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE CENTRE2(X)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION XMASS, YMASS, ZMASS, DIST, DISTMAX, X(3*NATOMS)
      INTEGER I, J1

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      IF (RIGID) THEN
         DO I=1,NATOMS/2
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=2*XMASS/NATOMS
         YMASS=2*YMASS/NATOMS
         ZMASS=2*ZMASS/NATOMS
C        PRINT*,'initial coordinates in centre:'
C        WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)
         DO I=1,NATOMS/2
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
      ELSE
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
         DO I=1,NATOMS
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(*,'(A,3F12.4)') 'centre of mass reset to the origin from ',XMASS,YMASS,ZMASS
C     PRINT*,'final coordinates in centre:'
C     WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)
C
C  Check that all the atoms are in the container. If not then rescale.
C  Must not do this - it could change the system from a minimum!
C
C     IF (RIGID) RETURN
C     DISTMAX=0.0D0
C     DO J1=1,NATOMS
C        DIST=X(3*J1-2)**2+X(3*J1-1)**2+X(3*J1)**2
C        IF (DIST.GT.DISTMAX) DISTMAX=DIST
C     ENDDO
C     IF (DISTMAX.GT.RADIUS) THEN
C        DISTMAX=DSQRT(DISTMAX/RADIUS)*0.99D0
C        DO J1=1,NATOMS
C           X(3*J1-2)=X(3*J1-2)*DISTMAX
C           X(3*J1-1)=X(3*J1-1)*DISTMAX
C           X(3*J1)  =X(3*J1)  *DISTMAX
C        ENDDO
C     ENDIF
      RETURN
      END
