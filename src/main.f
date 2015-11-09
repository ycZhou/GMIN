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
      PROGRAM GMIN
      use noa
      USE commons
      use qmodule
      use permu
      use f1com
      use modamber
      USE MODCHARMM
      IMPLICIT NONE

      INTEGER J1,J2, JP
      DOUBLE PRECISION, allocatable :: SCREENC(:)
      DOUBLE PRECISION POTEL

      COMMON /MYPOT/ POTEL

      CALL CPU_TIME(TSTART)
      call countatoms
      call modcommoninit
 
      allocate(fin(3*natoms))
      allocate(xicom(3*natoms),pcom(3*natoms))
      allocate(SCREENC(3*natoms))
      allocate(IATNUM(natoms), VT(natoms), ZSYM(natoms))
      VT(1:NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      IF (Q4T) CALL SHINIT

      CALL KEYWORD
      IF(CHRIGIDTRANST.OR.CHRIGIDROTT) CALL CHSETSEG
      IF(RMST) THEN
        ALLOCATE(RMSBEST(RMSSAVE,2),RMSCOOR(RMSSAVE,3*NATOMS))
        RMSBEST(1:RMSSAVE,1)=RMSLIMIT+RMSTOL
        RMSBEST(1:RMSSAVE,2)=0.D0
        RMSCOOR(1:RMSSAVE,1:3*NATOMS)=0.D0
        ALLOCATE(COORCOMP(1:3*NATOMS))
        OPEN(UNIT=1,FILE='compare',STATUS='OLD')
        READ(1,*) (COORCOMP(J1),J1=1,3*NATOMS)
        CLOSE(1)
      ENDIF
      allocate(FF(NSAVE),QMIN(NSAVE))
      allocate(QMINP(NSAVE,3*natoms))
      IF (GAUSST) THEN
         ALLOCATE(GAUSSKK(3*NATOMS,GMODES),GAUSSEE(GMODES),GKSMALL(3*NATOMS))
         CALL KEGEN ! initial setup 
         DO J1=1,GMODES
             PRINT *,J1,GAUSSEE(J1)
         ENDDO
      ENDIF
      QMINP(1:NSAVE,1:3*natoms)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*natoms,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memory
      VATO(1:natoms,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      ALLOCATE(ESAVE(NTAB,NPAR),XINSAVE(NTAB,NPAR))
      ALLOCATE(VEC(NVEC))
      IF (SYMMETRiZE.AND.(.NOT.CENT)) THEN
         PRINT '(A)','Probable input error - SYMMETRIZE true but CENT false'
         STOP
      ENDIF

      CALL IO1

      IF (SEEDT) THEN
         CALL GSEED
      ELSE
         IF ((.NOT.FIELDT).AND.CENT) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRE(COORDS,J1)
            ENDDO
         ENDIF
      ENDIF
      IF (SUPERSTEP) NSUPERSTEP=0
      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
      IF (NRUNS.GT.0) CALL MCRUNS(SCREENC)
C     CALL SYSTEM('rm ssdump ssave >& /dev/null')

      deallocate(fin)
      deallocate(xicom,pcom)
      IF (ALLOCATED(GAUSSKK)) DEALLOCATE(GAUSSKK,GAUSSEE)
      call modcommondeinit
      END
