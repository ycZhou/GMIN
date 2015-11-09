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
MODULE COMMONS
      use noa
      IMPLICIT NONE
      SAVE
     INTEGER :: NATOMS, NACCEPT, MAXIT, NTARGETS, NN, MM, NPAR, TIPID, NSYMREM=0, &
     &        NINTV, NSEED, NVEC, NS, NSSTOP, MAXIT2, NSAVE, &
     &        MCSTEPS(3), NRUNS, NVAR, MYPOWER, XMUL,  NTAB, NCOOP, NSYMINTERVAL, &
     &        NSUPER, NSACCEPT, NSUPERSTEP, GMODES, &
     &        NEVL, NEVS, NVECTORS, MUPDATE, LSYS, NHSMOVE, NRENORM, NRELAX, NHSRESTART, NRENSTUCK,&
     &        HBINS, NRBSITES, GATOM, NMSBSAVE, MAXSAVE, NSYMQMAX, EQUIL, NINTS, NumberOfTrajectories, &
     &        MCCycles, MCCyclesEquilibration, TargetWL, NTempPoints, SaveNth, DumpEveryNthQuench, NSpline, Hwindows, & 
     &        lhbins, sampledbins
 
     DOUBLE PRECISION RHO, GAMMA, SIG, SCEPS, SCC, TOLB, T12FAC, XMOVERENORM, RESIZE, QTSALLIS, &
     &                 CQMAX, RADIUS, BQMAX,  MAXBFGS, DECAYPARAM, SYMTOL1, SYMTOL2, SYMTOL3, SYMTOL4, SYMTOL5, &
     &                 ECONV, TOLD, TOLE, SYMREM(120,3,3), GMAX, CUTOFF, PCUT, EXPFAC, EXPD, &
     &                 BOXLX, BOXLY, BOXLZ, SUPSTEP, SQUEEZER, SQUEEZED, COOPCUT, STOCKMU, STOCKLAMBDA, &
     &                 TFAC(3), RMS, TEMPS, SACCRAT, CEIG, PNEWJUMP, EAMP, DISTFAC, &
     &                 APP, AMM, APM, XQP, XQM, ALPHAP, ALPHAM, ZSTAR, COMP, DGUESS, GUIDECUT, EFAC,& 
     &                 TRENORM, HISTMIN, HISTINT, HISTFAC, EPSSPHERE, FINALCUTOFF, &
     &                 HISTFACMUL, HPERCENT, AVOIDDIST, MAXERISE, TSTART, MATDIFF, STICKYSIG,&
     &                 MinimalTemperature, MaximalTemperature, SwapProb, hdistconstraint, &
     &                 RK_R, RK_THETA,ARMA,ARMB, ExtrapolationPercent, lnHarmFreq
      
     LOGICAL DEBUG, TARGET, MORSET, CUTT, SEEDT, CENT, TSALLIST, FREEZECORE, NEWJUMP, RENORM, CAPSID,& 
     &        OTPT, LJMFT, STRANDT, PAHT, SWT, MSTRANST, STOCKT, STICKYT, BLNT, &
     &        MSORIGT, SQUEEZET, PERIODIC, SCT, RESIZET, TIP, RIGID, Q4T, &
     &        SORTT, HIT, SAVEQ, PARALLELT, FIXD, RKMIN, BSMIN, PERMUTE, HIST, HISTRESTART, &
     &        SYMMETRIZE, DUMPT, NEON, ARGON, P46, NORESET, TABOOT, EVSTEPT, PACHECO, DL_POLY,&
     &        STAR, PLUS, TWOPLUS, GROUND, DIPOLE, DFTBT, SW, SUPERSTEP, EAMLJT, PBGLUET,&
     &        EAMALT, ALGLUET, MGGLUET, GUPTAT, DECAY, COOP, FIXBIN, GAUSST, &
     &        FRAUSIT, ANGST, SELFT, STEPOUT, WENZEL, THRESHOLDT, THOMSONT, &
     &        PROJ, RGCL2, TOSI, WELCH, AXTELL, AMBER, FIXIMAGE, BINARY, SHIFTCUT, ARNO, TUNNELT, TWOD,& 
     &        BLJCLUSTER, COMPRESST, FIX, FIXT, BFGS, LBFGST, DBRENTT, DZTEST, FNI, FAL, CPMD, TNT, ZETT1, &
     &        ZETT2, RESTART, CONJG, NEWRESTART, AVOID, NATBT, DIFFRACTT, CHRMMT, INTMINT, LB2T, PTMC, BINSTRUCTURES, &
     &        TETHER, HISTSMOOTH, VISITPROP, ARMT, FixedEndMoveT
    
      DOUBLE PRECISION :: DZP1, DZP2, DZP3, DZP4, DZP5, DZP6, DZP7
      LOGICAL :: FIELDT, OHT, IHT, TDT, D5HT
      DOUBLE PRECISION :: FOH, FIH, FTD, FD5H

      CHARACTER(LEN=80) :: SYS
      CHARACTER(LEN=1), ALLOCATABLE :: BEADLETTER(:), BLNSSTRUCT(:)

      DOUBLE PRECISION :: HESS(1,1)
!     DOUBLE PRECISION HESS(3*MXATMS,3*MXATMS)

!   allocatables

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NT, NQ, BLOCK, JUMPINT, JUMPTO  ! dimension will be NPAR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: HISTVALS, LHISTVALS  ! dimension declared HBINS in keyword

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: TEMP, STEP, OSTEP, ASTEP, ACCRAT, EPREV ! dimension will be NPAR 
      DOUBLE PRECISION, ALLOCATABLE :: TARGETS(:) ! allocated NTARGETS in keyword.f
      DOUBLE PRECISION, ALLOCATABLE :: ESAVE(:,:) ! dimension will be NTAB,NPAR
      DOUBLE PRECISION, ALLOCATABLE :: XINSAVE(:,:) ! dimension will be NTAB,NPAR
      DOUBLE PRECISION, ALLOCATABLE :: HDIST(:), HWEIGHT(:) ! dimension declared HISTBINMAX in keyword
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:) ! dimension would be NVEC - not used

      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FIXBOTH, FIXTEMP, FIXSTEP, JUMPMOVE, JDUMP, TMOVE, OMOVE ! dimension will be NPAR
      LOGICAL, ALLOCATABLE ::  IGNOREBIN(:) ! will be HBINS

      INTEGER, ALLOCATABLE, DIMENSION(:) ::  IATNUM     
      INTEGER,ALLOCATABLE :: ANV(:,:,:)
      DOUBLE PRECISION ,ALLOCATABLE :: LJREP_BLN(:,:), LJATT_BLN(:,:), A_BLN(:), B_BLN(:), C_BLN(:), D_BLN(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: TCOORDS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS,COORDSO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  VAT,  VATO    
      DOUBLE PRECISION, ALLOCATABLE :: SITE(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MSBCOORDS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MSBE(:)
      DOUBLE PRECISION, ALLOCATABLE :: GAUSSKK(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: GAUSSEE(:), GKSMALL(:)
      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: ZSYM   
      REAL, ALLOCATABLE :: xvEX(:)  !  single precision!!!!

      contains

      subroutine modcommoninit
         implicit none
         
         natoms = Number_of_Atoms
         allocate( ANV(NATOMS,NATOMS,3))         
      end subroutine modcommoninit

      subroutine modcommondeinit
         implicit none
         
         deallocate(anv)
      end subroutine modcommondeinit
END MODULE COMMONS
