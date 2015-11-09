!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MINPERMDIST(COORDSA,COORDSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2)
IMPLICIT NONE

INTEGER NATOMS, NPERM
INTEGER J3, INVERT, NORBIT1, NORBIT2, PERM(NATOMS), NCHOOSE2
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD
LOGICAL DEBUG, PERIODIC,TWOD

! OPEN(UNIT=10,FILE='minpermdist.xyz',STATUS='UNKNOWN')
! WRITE(10,'(I6)') NATOMS
! WRITE(10,'(A)') 'coordsa A'
! DO J3=1,NATOMS
!     WRITE(10,'(A2,2X,3F20.10)') 'LA',COORDSA(3*(J3-1)+1),COORDSA(3*(J3-1)+2),COORDSA(3*(J3-1)+3)
! ENDDO
! WRITE(10,'(I6)') NATOMS
! WRITE(10,'(A)') 'coordsb'
! DO J3=1,NATOMS
!     WRITE(10,'(A2,2X,3F20.10)') 'LA',COORDSB(3*(J3-1)+1:3*(J3-1)+3)
! ENDDO

! standard orientation - so permutational isomers should be found in one go
NCHOOSE2=1

20 CALL MYORIENT(COORDSB,DUMMYB,NORBIT1,1,NORBIT2,NCHOOSE2,NATOMS,DEBUG)
! WRITE(10,'(I6)') NATOMS
! WRITE(10,'(A)') 'dummyb'
! DO J3=1,NATOMS
!     WRITE(10,'(A2,2X,3F20.10)') 'LA',DUMMYB(3*(J3-1)+1:3*(J3-1)+3)
! ENDDO

DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)

INVERT=1
30 DUMMY(1:3*NATOMS)=INVERT*DUMMYA(1:3*NATOMS)
CALL MYORIENT(DUMMY,DUMMYA,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG)
! WRITE(10,'(I6)') NATOMS
! WRITE(10,'(A)') 'dummya'
! DO J3=1,NATOMS
!     WRITE(10,'(A2,2X,3F20.10)') 'LA',DUMMYA(3*(J3-1)+1:3*(J3-1)+3)
! ENDDO
!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance in this case.
!
10 CALL MINPERM(NATOMS, DUMMYB, DUMMYA, BOXLX, BOXLY, BOXLZ, PERIODIC, PERM, DISTANCE, DIST2, WORSTRAD)
! 10 CALL BIPARTITE(NATOMS, DUMMYB, DUMMYA, PERM, DISTANCE, DIST2, WORSTRAD)
! PRINT '(A,2G20.10)','DISTANCE,DIST2=',DISTANCE,DIST2
DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(PERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(PERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(PERM(J3)-1)+3)
   IF (J3.NE.PERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') 'permute atoms ',J3,' and ',PERM(J3)
      NPERM=NPERM+1
   ENDIF
ENDDO
IF (DEBUG) WRITE(*,'(A,I6,A,G20.10)') 'distance after permuting ',NPERM,' pairs of atoms=',SQRT(DISTANCE)
! WRITE(10,'(I6)') NATOMS
! WRITE(10,'(A,F12.5)') 'dummya after permuting; dist=',DISTANCE
! DO J3=1,NATOMS
!     WRITE(10,'(A2,2X,3F20.10)') 'LA',DUMMYA(3*(J3-1)+1:3*(J3-1)+3)
! ENDDO
!
!  Distance minimisations with respect to Euler angles and centre-of-mass.
!  Coordinates in DUMMYA are reset by mind (second argument).
!  
IF (NPERM.NE.0) THEN 
   CALL MINDGMIN(DUMMYB,DUMMYA,NATOMS,DISTANCE,PERIODIC,TWOD)
   IF (DEBUG) WRITE(*,'(A,G20.10)')  'distance after mind=                           ',DISTANCE
   IF (DISTANCE.LT.0.1D0) THEN
!     CLOSE(10)
      RETURN
   ENDIF
   GOTO 10
ELSE
   IF (INVERT.EQ.1) THEN
!
!  Now try the enantiomer if necessary.
!
      IF (DEBUG) PRINT '(A)','inverting geometry for comparison with target'
      INVERT=-1
      GOTO 30
   ENDIF
ENDIF
IF (NCHOOSE2.LT.NORBIT2) THEN
   NCHOOSE2=NCHOOSE2+1
   IF (DEBUG) PRINT '(2(A,I5))','trying NCHOOSE2=',NCHOOSE2,' for NORBIT2=',NORBIT2
   GOTO 20
ENDIF

! CLOSE(10)

END SUBROUTINE MINPERMDIST
