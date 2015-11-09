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
C  Energy and gradient for Thomson problem in theta, phi coordinates.
C
      SUBROUTINE THOMSON(X,V,ETHOMSON,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(*), DIST, V(*), ETHOMSON, DUMMY, CT1, ST1, CT2, ST2, CPDIFF, SPDIFF
      DOUBLE PRECISION COST(NATOMS), SINT(NATOMS), COSP(NATOMS), SINP(NATOMS)
      DOUBLE PRECISION, PARAMETER :: SR2=1.4142135623730950488D0

      DO J1=1,NATOMS
         J3=2*J1
         COST(J1)=COS(X(J3-1))
         SINT(J1)=SIN(X(J3-1))
         COSP(J1)=COS(X(J3))
         SINP(J1)=SIN(X(J3))
      ENDDO

      ETHOMSON=0.0D0
      VT(1:NATOMS)=0.0D0
      V(1:2*NATOMS)=0.0D0
      DO J1=1,NATOMS
         J3=2*J1
         CT1=COST(J1)
         ST1=SINT(J1)
         DO J2=J1+1,NATOMS
            J4=2*J2
            CT2=COST(J2)
            ST2=SINT(J2)
            CPDIFF=COSP(J1)*COSP(J2)+SINP(J1)*SINP(J2)
            SPDIFF=COSP(J2)*SINP(J1)-SINP(J2)*COSP(J1)
C           DIST=1.0D0/(SR2*SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2))
            DIST=1.0D0/SQRT(1.0D0-CT1*CT2-CPDIFF*ST1*ST2)
            ETHOMSON=ETHOMSON+DIST
            DIST=DIST**3
            DUMMY=SPDIFF*ST1*ST2*DIST
            V(J3-1)=V(J3-1)+(CPDIFF*CT1*ST2-CT2*ST1)*DIST
            V(J3)=V(J3)    -DUMMY
            V(J4-1)=V(J4-1)+(CPDIFF*CT2*ST1-CT1*ST2)*DIST
            V(J4)=V(J4)    +DUMMY
         ENDDO
      ENDDO
      ETHOMSON=ETHOMSON/SR2
      DO J1=1,2*NATOMS
         V(J1)=V(J1)/SR2
      ENDDO

      RETURN
      END
