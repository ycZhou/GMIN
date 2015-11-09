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
C     This subprogram performs a sort on the input data and
C     arranges it from smallest to biggest. The exchange-sort
C     algorithm is used.
C
      SUBROUTINE GSORT2(N,J3,A,B,NATOMS,F,NSAVE)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NATOMS, F(J3), NTEMP, NSAVE
      DOUBLE PRECISION TEMP, A(J3), B(NSAVE,3*NATOMS), C
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=F(L)
         F(L)=F(J1)
         F(J1)=NTEMP
         DO J2=1,3*NATOMS
            C=B(L,J2)
            B(L,J2)=B(J1,J2)
            B(J1,J2)=C
         ENDDO
20    CONTINUE
      RETURN
      END
