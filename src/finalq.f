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
      SUBROUTINE FINALQ
      USE commons
      use qmodule
      IMPLICIT NONE


      INTEGER J1, J2, ITERATIONS, BRUN,QDONE
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMS),
     1                  X(3*NATOMS), ENERGY,
     2                 GRAD(3*NATOMS), TIME

      COMMON /MYPOT/ POTEL
 
C
C  Make sure the lowest minima are tightly converged and then sort
C  them just to be on the safe side.
C
      IF (CUTT) CUTOFF=FINALCUTOFF
      SAVEQ=.FALSE.
      NQ(1)=0
      IF (FIELDT) FIELDT=.FALSE.
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      IF (SQUEEZET) SQUEEZET=.FALSE.
      MAXIT=MAXIT2
      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D6) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=QMINP(J1,J2)
            ENDDO
            NQ(1)=NQ(1)+1
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(*,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ(1),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO
         ENDIF
      ENDDO

      CALL GSORT2(NSAVE,NSAVE,QMIN,QMINP,NATOMS,FF,NSAVE)
C
C  Optionally sort the atoms from most bound to least bound according to VAT.
C
      IF (SORTT) THEN
         DO J1=1,NSAVE
            IF (QMIN(J1).LT.0.0D0) THEN
               DO J2=1,3*NATOMS
                  COORDS(J2,1)=QMINP(J1,J2)
                  X(J2)=QMINP(J1,J2)
               ENDDO
               CALL POTENTIAL(X,GRAD,ENERGY,.FALSE.,.FALSE.)
               DO J2=1,NATOMS
                  VAT(J2,1)=VT(J2)
               ENDDO
               CALL SORT3(NATOMS,NATOMS,VAT,COORDS,1,NPAR)
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=COORDS(J2,1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

C     IF (DEBUG) THEN
         IF (TABOOT) THEN
            IF (NPAR.GT.1) THEN
               PRINT*,'Taboo lists:'
               DO J1=1,NPAR
                  PRINT*,'Parallel run ',J1
                  WRITE(*,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
               ENDDO
            ELSE
               PRINT*,'Taboo list:'
               WRITE(*,'(6F15.7)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
C     ENDIF

      RETURN
      END
