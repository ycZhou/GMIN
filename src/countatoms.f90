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
MODULE NOA
      implicit none
      save
      integer :: Number_of_Atoms

      contains

      subroutine countatoms
      implicit none
      integer :: eof
      LOGICAL :: YESNO, YESNOA, YESNOC
      CHARACTER(LEN=10)  check
      CHARACTER(LEN=80) myline

      YESNO=.FALSE.
      YESNOA=.FALSE.
!
!  If the current working directory contains more than one of these files
!  then the precedence is coords, then input.crd, then coords.amber
!  OPTIM does this a bit better by calling getparams first to see if
!  we are actually doing AMBER or CHARMM. 
!
      INQUIRE(FILE='coords',EXIST=YESNO)
      INQUIRE(FILE='coords.amber',EXIST=YESNOA)
      INQUIRE(FILE='input.crd',EXIST=YESNOC)

      Number_of_Atoms=0

      IF (YESNO) THEN
         OPEN(UNIT=7,FILE='coords',STATUS='OLD')
         PRINT '(A)','reading coordinates from file coords'
         do
            read(7,*,iostat=eof)
            if (eof==0) then
               Number_of_Atoms = Number_of_Atoms + 1
            else
               exit
            endif
         enddo
      ELSEIF (YESNOC) THEN
         PRINT '(A)','reading coordinates from file input.crd'
         OPEN(UNIT=7,FILE='input.crd',STATUS='OLD')
         do
           read(7,*) myline
           if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
              cycle
           else
              read(myline,*) Number_of_Atoms
              exit
           endif
         enddo

! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.

         CALL GETMAXAIM
      ELSEIF (YESNOA) THEN
         OPEN(UNIT=7,FILE='coords.amber',STATUS='OLD')
         PRINT '(A)','reading coordinates from file coords.amber'
         do
            read(7,'(A3)',iostat=eof) check
            if (eof.LT.0) then
               PRINT *,'End of file before all information specified'
               STOP
            ENDIF
            IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') THEN
               CLOSE(7)
               EXIT
            ENDIF
            Number_of_Atoms = Number_of_Atoms + 1
         enddo
      ELSE
         PRINT '(A)','ERROR - no coords, input.crd or coords.amber file'
         STOP
      ENDIF

      CLOSE(7)

      print *, Number_of_Atoms, ' atoms in the system'
      
      end subroutine countatoms

END MODULE NOA
