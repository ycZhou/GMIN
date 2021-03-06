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
SUBROUTINE ENUMERATE(NFLOAT,NEWORB,NEWORBSIZE,HARDMAXIMUM,OCCS,NPOSS,DEBUG)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NFLOAT, NEWORB, HARDMAXIMUM
INTEGER, INTENT(IN) :: NEWORBSIZE(NEWORB)
INTEGER COUNTS(0:NFLOAT)
INTEGER J, I, MAXCOUNT, MAXDEGEN, J1, NDUMMY, NOK
INTEGER, ALLOCATABLE :: OCCUPATIONS(:,:,:), SAVEOCC(:,:,:)
INTEGER, INTENT(OUT) :: OCCS(HARDMAXIMUM,NEWORB), NPOSS
LOGICAL DEBUG

MAXDEGEN=MIN(HARDMAXIMUM,100)
! ALLOCATE(OCCUPATIONS(0:NFLOAT,MAXDEGEN,NEWORB))
ALLOCATE(OCCUPATIONS(0:NFLOAT,HARDMAXIMUM,NEWORB))
COUNTS(0)=1; COUNTS(1:NFLOAT)=0
MAXCOUNT=0
OCCUPATIONS(0:NFLOAT,1:MAXDEGEN,1:NEWORB)=0

IF (DEBUG) PRINT '(A,3I6)','in enumerate, NFLOAT, HARDMAXIMUM, NEWORB=',NFLOAT, HARDMAXIMUM, NEWORB

DO J=1,NEWORB
   DO I=NEWORBSIZE(J),MIN(NFLOAT,NEWORBSIZE(J)+MAXCOUNT)
      NOK=0
      DO J1=1,COUNTS(I-NEWORBSIZE(J))
         IF (OCCUPATIONS(I-NEWORBSIZE(J),J1,J).EQ.0) THEN
            NOK=NOK+1

!           IF ((I.EQ.NFLOAT).AND.(COUNTS(I)+NOK.GT.HARDMAXIMUM)) THEN
!              OCCS(1:COUNTS(I),1:NEWORB)=OCCUPATIONS(NFLOAT,1:COUNTS(I),1:NEWORB)
!              NPOSS=COUNTS(I)
!              RETURN
!           ENDIF
               
!           IF (COUNTS(I)+NOK.GT.MAXDEGEN) THEN
!              ALLOCATE(SAVEOCC(0:NFLOAT,MAXDEGEN,NEWORB))
!              SAVEOCC(0:NFLOAT,1:MAXDEGEN,1:NEWORB)=OCCUPATIONS(0:NFLOAT,1:MAXDEGEN,1:NEWORB)
!              DEALLOCATE(OCCUPATIONS)
!              MAXDEGEN=MAXDEGEN*2
!              ALLOCATE(OCCUPATIONS(0:NFLOAT,MAXDEGEN,NEWORB))
!              OCCUPATIONS(0:NFLOAT,1:MAXDEGEN/2,1:NEWORB)=SAVEOCC(0:NFLOAT,1:MAXDEGEN/2,1:NEWORB)
!              OCCUPATIONS(0:NFLOAT,MAXDEGEN/2+1:MAXDEGEN,1:NEWORB)=0
!              DEALLOCATE(SAVEOCC)
!              PRINT '(A,I10)','allocation size increased to ',MAXDEGEN
!           ENDIF

            IF (COUNTS(I)+NOK.LE.HARDMAXIMUM) THEN
               OCCUPATIONS(I,COUNTS(I)+NOK,1:NEWORB)=OCCUPATIONS(I-NEWORBSIZE(J),J1,1:NEWORB)
               OCCUPATIONS(I,COUNTS(I)+NOK,J)=OCCUPATIONS(I-NEWORBSIZE(J),J1,J)+1
            ENDIF
         ENDIF
      ENDDO
      COUNTS(I)=MIN(COUNTS(I)+NOK,HARDMAXIMUM)
   ENDDO
   MAXCOUNT=MAXCOUNT+NEWORBSIZE(J)
ENDDO

IF (DEBUG) PRINT '(A)','final counts:'
IF (DEBUG) PRINT '(2I16)',(J,COUNTS(J),J=0,NFLOAT)

OCCS(1:COUNTS(NFLOAT),1:NEWORB)=OCCUPATIONS(NFLOAT,1:COUNTS(NFLOAT),1:NEWORB)
NPOSS=COUNTS(NFLOAT)
IF (DEBUG) PRINT '(A,I5)','number of possibilities=',NPOSS

RETURN

! DO J=0,NFLOAT
!    IF (COUNTS(J).GT.0) THEN
!       PRINT '(A,I5,A)','occupations for ',J,' atoms:'
!       DO I=1,COUNTS(J)
!          NDUMMY=0
!          DO J1=1,NEWORB
!             NDUMMY=NDUMMY+OCCUPATIONS(J,I,J1)*NEWORBSIZE(J1)
!          ENDDO
!          PRINT '(20I4)',OCCUPATIONS(J,I,1:NEWORB)
!          PRINT '(A,I5)','checksum=',NDUMMY
!       ENDDO
!    ENDIF
! ENDDO

RETURN
END
