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
C  Program to convert capsid CofM and DV coordinates to pentagons.
C
      PROGRAM CAPSID6to15
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X1, Y1, Z1, 
     2                 M1, L1, N1, ALPHA1, RAD, CA1, SA1, C3A1,
     3                 NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

      RAD=1.1D0
      NUM1=-(1.0D0+SQRT(5.0D0))*RAD/4.0D0
      NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)*RAD/2.0D0
      NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)*RAD/2.0D0
      NUM4=(SQRT(5.0D0)-1.0D0)*RAD/4.0D0
      NUM5=-(1.0D0+SQRT(5.0D0))*RAD

10    READ(*,*,END=666) X1, Y1, Z1, L1, M1, N1
      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=-ALPHA1/2+ALPHA1**3/24
         C3A1=-0.5D0+ALPHA1**2/24.0D0
         SA1=1.0D0-ALPHA1**2/6
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         SA1=SIN(ALPHA1)/ALPHA1
      ENDIF

      WRITE(*,'(3F20.10)')
     @     x1 + rad*CA1 - l12*rad*C3A1,
     @     y1 - l1*m1*rad*C3A1 - n1*rad*SA1,
     @     z1 - l1*n1*rad*C3A1 + m1*rad*SA1,
     @     x1 + NUM4*CA1 - l1*(m1*NUM3 + l1*NUM4)*C3A1 + n1*NUM3*SA1,
     @     y1 + NUM3*CA1 - m1*(m1*NUM3 + l1*NUM4)*C3A1 - n1*NUM4*SA1,
     @     z1 - n1*(m1*NUM3 + l1*NUM4)*C3A1 + (-(l1*NUM3) + m1*NUM4)*SA1,
     @     x1 + NUM1*CA1 - l1*(l1*NUM1 + m1*NUM2)*C3A1 + n1*NUM2*SA1,
     @     y1 + NUM2*CA1 - m1*(l1*NUM1 + m1*NUM2)*C3A1 - (n1*NUM5*SA1)/4.,
     @     z1 - n1*(l1*NUM1 + m1*NUM2)*C3A1 + (m1*NUM1 - l1*NUM2)*SA1,
     @     x1 + NUM1*CA1 - l1*(l1*NUM1 - m1*NUM2)*C3A1 - n1*NUM2*SA1,
     @     y1 - NUM2*CA1 - m1*(l1*NUM1 - m1*NUM2)*C3A1 - (n1*NUM5*SA1)/4.,
     @     z1 - n1*(l1*NUM1 - m1*NUM2)*C3A1 + (m1*NUM1 + l1*NUM2)*SA1,
     @     x1 + NUM4*CA1 - l1*(-(m1*NUM3) + l1*NUM4)*C3A1 - n1*NUM3*SA1,
     @     y1 - NUM3*CA1 - m1*(-(m1*NUM3) + l1*NUM4)*C3A1 - n1*NUM4*SA1,
     @     z1 - n1*(-(m1*NUM3) + l1*NUM4)*C3A1 + (l1*NUM3 + m1*NUM4)*SA1

      GOTO 10

666   STOP
      END
