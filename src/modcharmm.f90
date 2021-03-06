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
MODULE modcharmm
   IMPLICIT NONE
   SAVE

      INTEGER, DIMENSION(:), ALLOCATABLE  :: NPHIPSI,NOMEGAC,NSIDECHAIN,NCHIRAL  ! CHDIHE
      INTEGER, DIMENSION(:), ALLOCATABLE :: PHIPSI,OMEGAC,SIDECHAIN,CHIRAL  ! CHDIHE
      INTEGER :: CHARMMTYPE, CHSEEDI, CHSEEDJ, CHSEEDK    ! CHARMMWORD
      INTEGER :: CHFREQ, FTRANS, FROT, RMSSAVE, NEWCONFST
      LOGICAL :: TOMEGAC,TSIDECHAIN     ! CHARMMTWIST
      LOGICAL :: DAESTAT, NOPHIPSIT, OMEGAT, NEWCONFT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: REFCOORD,COORCOMP   ! CHREF
      DOUBLE PRECISION :: CHPMIN, CHPMAX, CHNMIN, CHNMAX
      DOUBLE PRECISION :: PTRANS, PROT, TRANSMAX, ROTMAX, RMSLIMIT, RMSTOL, NCWALL
      INTEGER ISEED
      LOGICAL :: GCHARMMFAIL,CHRIGIDTRANST, CHRIGIDROTT, RMST, SELECTT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSEGATOMS     ! dimension will be NSEG
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEGMENTT
      DOUBLE PRECISION :: PIVOTP, SIDECHAINP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: RMSBEST,RMSCOOR 

END MODULE
