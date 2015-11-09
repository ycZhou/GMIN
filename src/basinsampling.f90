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
! 
! The `basin-sampling' algorithm combines `basin-hopping' and Wang-Landau
! sampling technique to study the thermodynamics of the transformed PES.
! It provides a direct temperature-independent estimate of the density of states
! along with thermodynamic properties such as the free energy and entropy via ensemble averages using
! samples of local minima rather than instantaneous configurations. (Tetyana Bogdan)
! T.V. Bogdan, D.J. Wales and F. Calvo, J. Chem. Phys., 124, 044102 (2006).
!                                            
!---======================================---
      subroutine BasinSampling

      use modcharmm
      use Commons, only: Natoms,  coords, histrestart, histmin, histint, hbins, hist, histfac, TargetWL, &
                         & histfacmul, debug, Equil , periodic, twod, binstructures, SaveNth, DumpEveryNthQuench, &
                         & FixedEndMoveT, chrmmt 
      implicit none
    

      integer nsteps, i, j, Visits(Hbins), VisitsTotal(Hbins), lbfgs_iterations, BinIndexold, BinIndexnew, NQuenches, &
              & NQuenchesSinceFlat, NQuenchesSuccess, nWL, NQuenchesSinceLastUpdate, MinimaNumber(Hbins), ndummy, Converged, &
              & FrozenVisits(Hbins)
      real(8) potel, dummy, CurrentPointEnergy, CurrentPointCoordinates(3*Natoms, 1), &
              & Distanceold, lnModFac, PerturbedCoordinates(3*Natoms, 1), PerturbedCoordinatesSave(3*Natoms, 1), lnRatio, &
              & harvest, BinLabel(HBins), oldenergy, Norm, MinDistance(Hbins), MinDistanceold, lnharvest, QWeight(hbins,2), &
              & CurrentQ, Q4orderparam, knownenergies(61836)

      !real(16) lnWeight(Hbins), Distance(Hbins) 
      real(8) lnWeight(Hbins), Distance(Hbins) 
         
      logical yesno, Flat, evap, evapreject, AcceptMove

      common /mypot/ potel

      common /ev/ evap, evapreject

! If the BS run is requested to be restarted the following files will be read in:
      
      if (histrestart) then
         open(unit=421, file='lnWeight.his', status='old')
         do i=1, Hbins
            read(421, '(2G20.10)') dummy, lnWeight(i)
         enddo
         print *, 'Following ln(Weight) histogram was read:'
         print *, lnWeight
         close(421)
         open(unit=421, file='Distance.his', status='old')
         do i=1, Hbins
            read(421, '(2G20.10)') dummy, Distance(i)
         enddo
         print *, 'Following Distance histogram was read:'
         print *, Distance
         close(421)
         open(unit=421, file='MinDistance.his', status='old')
         do i=1, Hbins
            read(421, '(2G20.10)') dummy, MinDistance(i)
         enddo
         print *, 'Following Minimized Distance histogram was read:'
         print *, MinDistance
         close(421)
      else
! note that weight accumulation is is done in logarithms
         lnWeight=0.0d0 
         QWeight=0.0d0
         Distance=0.0d0
         MinDistance=0.0d0
      endif

      if (histrestart) then
         open(unit=422, file='Visits.his', status='old')
         do i=1, Hbins
            read(422, '(2G20.10)') dummy, Visits(i)
         enddo
         print *, 'Following Visits histogram was read:'
         print *, Visits
         close(422)
         open(unit=422, file='VisitsTotal.his', status='old')
         do i=1, Hbins
            read(422, '(2G20.10)') dummy, VisitsTotal(i)
         enddo
         print *, 'Following VisitsTotal histogram was read:'
         print *, VisitsTotal
         close(422)
         open(unit=422, file='lnModfac.restart', status='old')
         read(422, '(G20.10)') lnModFac
         print *, 'Restarting from modification factor', lnModFac
         close(422)
         open(unit=422, file='nWL.restart', status='old')
         read(422, '(G20.10)') nWL
         print *, 'Number of Wang Landau iterations already completed: ', nWL
         close(422)
      else
         Visits=0
         VisitsTotal=0
         MinimaNumber=0
         lnModFac=log(HistFac) ! refer to the Table
         print *, 'lnModfac, exp', lnModFac, exp(lnModFac)
         nWL=0
      endif

      do i=1, Hbins
         BinLabel(i)=HistMin + HistInt*(i-0.5)
      enddo


! for VISITPROP convergence scenario 
! if one wants to use previously accumulated number of visits, set Equil to -1 in data and copy 
! old VisitsTotal.his to Frozen.his
!
      if (Equil.eq.-1) then
         open(unit=422, file='Frozen.his', status='old')
         do i=1, Hbins
            read(422, '(2G20.10)') dummy, FrozenVisits(i)
         enddo
       endif

      
      CurrentPointCoordinates=coords
      PerturbedCoordinatesSave=coords
      print *, 'Calculating initial energy'
      call quench(.false.,1,lbfgs_iterations,dummy,ndummy,Converged,CurrentPointCoordinates)
      CurrentPointEnergy=potel
      print *, 'Initial energy', CurrentPointEnergy


      BinIndexold=Energy2Index(CurrentPointEnergy)

      if ((BinIndexold < 1 ).or.(BinIndexold > HBins)) then
         print *, 'Starting geometry is outside requested range. Exiting.'
         return
      endif 
      Distanceold=0.0d0       

      NQuenches=0
      NQuenchesSinceFlat=0
      NQuenchesSuccess=0
      Converged=0
      NQuenchesSinceLastUpdate=0

! Repeat for the requested number of Wang-Landau iterations. 
! The iteration is complete when the visits histogram satisfies the flatness criterion.

      do
         if ( nWL==TargetWL ) exit
! jmc do something different here if charmm
         IF (CHRMMT) THEN

! Fixed end move scheme
             IF (FixedEndMoveT) THEN
!!!                 CALL FixedEndMove(1)
!!!                 PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
!!!             ELSE
! step in internal coordinates
                ! first make sure the internal coordinate table corresponds to the starting Cartesians
                ! (important if the previous quench failed...)
                coords(:3*NATOMS,1) = PerturbedCoordinatesSave(:3*NATOMS,1)
                CALL FILLICT(1)
                IF(CHRIGIDTRANST) CALL MKRIGIDTRANS(1)
                IF(CHRIGIDROTT) CALL MKRIGIDROT(1)
                CALL TAKESTEPCH(1)
                ! new geometry is now in coords
                PerturbedCoordinates(:3*NATOMS,1) = coords(:3*NATOMS,1)
             ELSE ! stealing FEM keyword for basin-hopping gmin IC moves...
! try taking step in Cartesians; need to preserve detailed balance and explore all of configuration space.
                PerturbedCoordinates = PerturbGeometry(PerturbedCoordinatesSave)
                coords(1:3*NATOMS,1) = PerturbedCoordinates(1:3*NATOMS,1)
                CALL FILLICT(1)
             ENDIF
         ELSE
         PerturbedCoordinates=PerturbGeometry(PerturbedCoordinatesSave) 
         coords=PerturbedCoordinates
         ENDIF

         call quench(.false.,1,lbfgs_iterations,dummy,ndummy,Converged,coords)
         if (evapreject) then
            if (debug) print *, 'Evaporation during minimisation'
            if (debug) print *, 'oldenergy, CurrentPointEnergy, potel=',oldenergy, CurrentPointEnergy, potel
            cycle
         endif
         oldenergy=CurrentPointEnergy
         CurrentPointEnergy=potel
         CurrentPointCoordinates=coords


         if ((Converged.eq.1).and.(.not.evapreject)) then
            NQuenchesSuccess=NQuenchesSuccess+1
            BinIndexnew=Energy2Index(CurrentPointEnergy)
            !call ORDERQ4(Natoms, CurrentPointCoordinates, Q4orderparam)
            !CurrentQ=Q4orderparam
            !print *, BinIndexnew, CurrentPointEnergy, CurrentQ, Q4orderparam


            if (BinIndexnew < 1 .or.BinIndexnew>Hbins) then
               if (debug) print *, 'Structure outside energy range. Rejecting move'
               AcceptMove=.false.
            else 
         	if ((binstructures).and.(mod(NQuenches, SaveNth).EQ.0)) then
             	   call SaveBinStructures(CurrentPointEnergy, CurrentPointCoordinates, BinIndexnew, MinimaNumber, .true.)
           	 endif
               call SaveBinStructures(CurrentPointEnergy, CurrentPointCoordinates, BinIndexnew, MinimaNumber, .false.)
               lnRatio=lnWeight(BinIndexold)-lnWeight(BinIndexnew)
               call random_number(harvest)
               lnharvest=log(harvest)
               if ((lnRatio>0.d0).or.(lnRatio>lnharvest)) then
                  AcceptMove=.true.
               else
                  AcceptMove=.false.
               endif
            endif

            if (debug) print *, 'Eold, Enew, Iold, Inew, Wold, Wnew, Converged , Ratio'
            if (debug) print *, oldenergy, CurrentPointEnergy, BinIndexold, BinIndexnew, lnWeight(BinIndexold) , &
                     & lnWeight(BinIndexnew), Converged, lnRatio

            if (AcceptMove) then 
               ! move accepted
               !print *, 'Moving from bin', BinIndexold, 'to bin', BinIndexnew

               Visits(BinIndexnew)=Visits(BinIndexnew)+1
               VisitsTotal(BinIndexnew)=VisitsTotal(BinIndexnew)+1
               lnWeight(BinIndexnew)=lnWeight(BinIndexnew)+lnModFac
               Distanceold=CalculatedDistance(CurrentPointCoordinates, PerturbedCoordinates)
               call mindGMIN(CurrentPointCoordinates,PerturbedCoordinates,natoms,MinDistanceold,periodic,twod)
               Distance(BinIndexnew)=Distance(BinIndexnew)+Distanceold
               MinDistance(BinIndexnew)=MinDistance(BinIndexnew)+MinDistanceold
               BinIndexold=BinIndexnew
               PerturbedCoordinatesSave=PerturbedCoordinates
            else
               ! move rejected
               !print *, 'Staying in bin', BinIndexold
               Visits(BinIndexold)=Visits(BinIndexold)+1
               VisitsTotal(BinIndexold)=VisitsTotal(BinIndexold)+1
               lnWeight(BinIndexold)=lnWeight(BinIndexold)+lnModFac
               Distance(BinIndexold)=Distance(BinIndexold)+Distanceold
               MinDistance(BinIndexold)=MinDistance(BinIndexold)+MinDistanceold
            endif
         else
         if (debug) print *, 'Optimisation unsuccessful'
         endif

         if ((mod(NQuenchesSuccess,DumpEveryNthQuench).eq.0).and.(NQuenchesSinceLastUpdate > Equil)) then  
            ! record statistics
             open(unit=421, file='lnWeight.his', status='unknown')
             do i=1, Hbins
                write(421, '(2G20.10)')  BinLabel(i), lnWeight(i)
             enddo
             close(421)
           
             open(unit=421, file='Distance.his', status='unknown')
             do i=1, Hbins
                write(421, '(2G20.10)')  BinLabel(i), Distance(i)
             enddo
             close(421)
           
             open(unit=421, file='MinDistance.his', status='unknown')
             do i=1, Hbins
                write(421, '(2G20.10)')  BinLabel(i), MinDistance(i)
             enddo
             close(421)
      
             open(unit=422, file='Visits.his', status='unknown')
             do i=1, Hbins
                write(422, '(2G20.10)')  BinLabel(i), Visits(i)
             enddo
             close(422)
      
             open(unit=422, file='VisitsTotal.his', status='unknown')
             do i=1, Hbins
                write(422, '(2G20.10)')  BinLabel(i), VisitsTotal(i)
             enddo
             close(422)
           
             open(unit=422, file='lnModfac.restart', status='unknown')
             write(422, '(G20.10)') lnModFac
             close(422)
           
             open(unit=422, file='nWL.restart', status='unknown')
             write(422, '(G20.10)') nWL
             close(422)
      
             open(unit=421, file='MinimaNumber.his', status='unknown')
             do i=1, Hbins
                write(421, '(2G20.10)')  BinLabel(i), MinimaNumber(i)
             enddo
             close(421)


             call Record_Stat(lnWeight,Distance, MinDistance, VisitsTotal, BinLabel)

             ! checking histogram for convergence

             Flat=CheckFlatness(Visits,lnModFac,nWL, FrozenVisits)
             if (Flat) then
                NQuenchesSinceLastUpdate=0
                print *, '--===Visits histogram satisfied flatness criterion===--'
                lnModFac=HistFacMul*lnModFac
                print *, 'lnModfac, exp', lnModFac, exp(lnModFac)
                print *, 'Updating modification factor to', lnModFac
                FrozenVisits=Visits
                Visits=0
                nWL=nWL+1
                
                open(unit=422, file='lnModfac.restart', status='unknown')
                write(422, '(G20.10)') lnModFac
                close(422)
              
                open(unit=422, file='nWL.restart', status='unknown')
                write(422, '(G20.10)') nWL
                close(422)
             endif
      

              print *, NQuenches, ' quenches completed'
         endif

      NQuenches=NQuenches+1
      NQuenchesSinceLastUpdate=NQuenchesSinceLastUpdate+1
      enddo

      contains

!---======================================---`
      function Energy2Index(CurrentPointEnergy)
      use Commons
      implicit none

      integer Energy2Index

      real(8) CurrentPointEnergy


      if (nint((CurrentPointEnergy-HistMin)/HistInt)<0) then
         Energy2Index=-1
      else
         Energy2Index=nint((CurrentPointEnergy-HistMin)/HistInt)+1
      endif

      end function Energy2Index


!---======================================---`
      real(8) function GetDisplacement(step)
           implicit none
           real(8),intent(in) :: step
           real(8) :: harvest
           call random_number(harvest)
           GetDisplacement=2.0d0*harvest*step-step
      end function GetDisplacement

!---======================================---`
      function PerturbGeometry(CurrentPointCoordinates)
           use Commons, only: Radius,Natoms,step
           implicit none

           real(8),dimension(3*Natoms,1) :: PerturbGeometry

           real(8),intent(in) :: CurrentPointCoordinates(3*Natoms,1)
           integer i, j

           PerturbGeometry=CurrentPointCoordinates

           do j=1, Natoms
                do
                     PerturbGeometry(3*(j-1)+1,1) = PerturbGeometry(3*(j-1)+1,1) + GetDisplacement(step(1))
                     PerturbGeometry(3*(j-1)+2,1) = PerturbGeometry(3*(j-1)+2,1) + GetDisplacement(step(1))
                     PerturbGeometry(3*(j-1)+3,1) = PerturbGeometry(3*(j-1)+3,1) + GetDisplacement(step(1))
                     if ( PerturbGeometry(3*(j-1)+1, 1)**2 &
                      & + PerturbGeometry(3*(j-1)+2, 1)**2 &
                      & + PerturbGeometry(3*(j-1)+3, 1)**2 < Radius ) then
                          exit
                     else
                          PerturbGeometry(3*(j-1)+1,1) = CurrentPointCoordinates(3*(j-1)+1,1)
                          PerturbGeometry(3*(j-1)+2,1) = CurrentPointCoordinates(3*(j-1)+2,1)
                          PerturbGeometry(3*(j-1)+3,1) = CurrentPointCoordinates(3*(j-1)+3,1)
                     endif
		enddo
           enddo
      end function PerturbGeometry

!---======================================---`

      function CalculatedDistance(CurrentPointCoordinates, PerturbedCoordinates)
      use Commons
      implicit none

      integer i
      real(8) CalculatedDistance, CurrentPointCoordinates(3*Natoms,1), PerturbedCoordinates(3*Natoms,1), &
              & x1(Natoms, 1), y1(Natoms, 1), z1(Natoms, 1), x2(Natoms, 1), y2(Natoms, 1), z2(Natoms, 1)

      
      do i=1,Natoms
         x1(i,1)=CurrentPointCoordinates(3*(i-1)+1,1)
         y1(i,1)=CurrentPointCoordinates(3*(i-1)+2,1)
         z1(i,1)=CurrentPointCoordinates(3*(i-1)+3,1)
         x2(i,1)=PerturbedCoordinates(3*(i-1)+1,1)
         y2(i,1)=PerturbedCoordinates(3*(i-1)+2,1)
         z2(i,1)=PerturbedCoordinates(3*(i-1)+3,1)
      enddo

      CalculatedDistance=dsqrt(sum ( (x1(:,1)-x2(:,1))**2 + (y1(:,1)-y2(:,1))**2 + (z1(:,1)-z2(:,1))**2 ) )

      if (debug) print *, 'Distance' , CalculatedDistance

      end function CalculatedDistance

!---======================================---`
      function CheckFlatness(Visits,lnModFac,nWL, FrozenVisits)
      use Commons
      implicit none

      integer Visits(Hbins), n, i, VisitsNonzero(Hbins),Hmin, Hmax, deltaHk, nWL, FrozenVisits(Hbins)
      logical CheckFlatness, FlatBin(Hbins)
      real(8) HDev, lnModFac
      character (len =256) filename
      character (len= 10)  istr


      CheckFlatness=.false.
      deltaHk=0


      Hmax=maxval(Visits)
      do i=1, Hbins
         if (Visits(i).ne.0) then
            VisitsNonzero(i)=Visits(i)
         else
            VisitsNonzero(i)=huge(hbins)
         endif
      enddo
      Hmin=minval(VisitsNonzero)
      if (debug) print *, 'Min and max', Hmin, Hmax

      HDev=(1.0d0*Hmax-Hmin)/(Hmax+Hmin)
    
      open(unit=421, file='HistFluctuations', status='unknown', position='append')
      write(421, '(G20.10)') HDev
      close(421)

      if (VisitProp) then
          !minimal number of visits should be proportional to 1/sqrt(ln(f))
 
          CheckFlatness=.false.
          FlatBin=.false.
 
          if ((Equil.ne.-1).and.(nWL.eq.0)) then
             if ((HDev<HPercent).and.(abs(1.0d0*Hmin-Hmax)>HPercent)) CheckFlatness=.true.
          else
             do i=1, Hbins
                if ((Visits(i) > 1.0d0/sqrt(lnModFac)).or.((Visits(i).eq.0).and.(FrozenVisits(i).eq.0))) FlatBin(i)=.True.
             enddo
             if (debug)  print *, 'Flatness by bin', FlatBin 
             if (All(FlatBin)) CheckFlatness=.true.
          endif
      else
         ! histogram flatness criterion: should be used for disconnected PES in the initial run
         ! where 1/lnf criterion will never be satisfied.
         if ((HDev<HPercent).and.(abs(1.0d0*Hmin-Hmax)>HPercent)) CheckFlatness=.true.
      endif

      do i=1, Hbins
         if (Visits(i).ne.0) then
            deltaHk=deltaHk+(Visits(i)-Hmin)
         else
            deltaHk=deltaHk+0
         endif
      enddo
      open(unit=421, file='HistSaturation', status='unknown', position='append')
      write(421, '(G20.10)') deltaHk
      close(421)

      write (istr, '(i10)') nWL
      filename="HistSaturation."//trim(adjustl(istr))
      open(unit=421,file=filename, status="unknown", form="formatted", position="append")
      write(421, '(G20.10)') deltaHk
      close(421)

      if (CheckFlatness) then
         open(unit=421, file='deltaHk.vs.lnModfac', status='unknown', position='append')
         write(421, '(2G20.10)') lnModFac,deltaHk
         close(421)
      endif

      end function CheckFlatness

!---======================================---`
      subroutine Record_Stat(lng_j,Distance, MinDistance, VisitsTotal, BinLabel)
      use Commons
      implicit none

      integer VisitsTotal(Hbins), i, t, fullbins

      real(8)  MinDistance(hbins), kB, BinLabel(Hbins),  &
              & a, b, ap, bp, avlng_j, totlng_j, minimal_d_j, &
              & deltaT, sumd_jnonmin, sumd_jmin, &
              & avd_jnonmin, avd_jmin, sumd2_jnonmin, sumd2_jmin, avd2_jnonmin, avd2_jmin, &
	      & std_devd_jnonmin, std_devd_jmin, std_devsmoothd_j, avsmoothd2_j, avsmoothd_j, &
              & totlnp_j, avlnp_j

      real(8) lnp_jmin(Hbins), p_jnonmin(Hbins), lnp_jnonmin(hbins), lng_j(Hbins), kappa, &
              & d_jnonmin(Hbins), Distance(Hbins), d_jmin(hbins), scaledlng_jnonmin(hbins),&
              & p_jnorm(Hbins), smoothd_jnonmin(hbins) , g_j(Hbins), g_jnorm(Hbins), p_j(Hbins), scaledlnp_jnonmin(hbins)


      if (periodic) then
         kappa=3*Natoms-3
      else
         kappa=3*Natoms-6
      endif

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            d_jnonmin(i)=Distance(i)/VisitsTotal(i)
         else
            d_jnonmin(i)=0.0d0
         endif
      enddo 

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            d_jmin(i)=MinDistance(i)/VisitsTotal(i)
         else
            d_jmin(i)=0.0d0
         endif
      enddo 

      minimal_d_j=1.0d100
      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
           if (d_jnonmin(i) < minimal_d_j) minimal_d_j=d_jnonmin(i)
         endif
      enddo

     ! calculate density of minima in the BS approximation: P_j=[G_j*(dm/d_j)**kappa]/Norm*deltaV 
     ! now we will actually be storing ln(p_j)

      fullbins=0
      totlng_j=0.d0
      avlng_j=0.d0

      do i=1, Hbins
         totlng_j=totlng_j+lng_j(i)
         fullbins=fullbins+1
      enddo

      avlng_j=totlng_j/fullbins

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            scaledlng_jnonmin(i)=lng_j(i)-avlng_j
         else
            scaledlng_jnonmin(i)=0.d0
         endif
      enddo
 
      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            g_j(i)=exp(scaledlng_jnonmin(i))
         else
            g_j(i)=0.d0
         endif
      enddo

	do i=1,Hbins
	    if (VisitsTotal(i).ne.0) then
	       g_jnorm(i)=g_j(i)/(sum(g_j)*HistInt)
	    else
	       g_jnorm(i)=0.d0
            endif
        enddo

         
      fullbins=0
      totlnp_j=0.d0
      avlnp_j=0.d0

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            lnp_jnonmin(i)=lng_j(i)-(kappa+3)*log(d_jnonmin(i))
            totlnp_j=totlnp_j+lnp_jnonmin(i)
            fullbins=fullbins+1 
         else 
            lnp_jnonmin(i)=-Huge(lng_j(i))
         endif
      enddo

      avlnp_j=totlnp_j/fullbins

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            scaledlnp_jnonmin(i)=lnp_jnonmin(i)-avlnp_j
         else
            scaledlnp_jnonmin(i)=0.d0
         endif
      enddo

      do i=1, Hbins
         if (VisitsTotal(i).ne.0) then
            p_jnonmin(i)=exp(scaledlnp_jnonmin(i))
         else
            p_jnonmin(i)=0.d0
         endif
      enddo

        do i=1,Hbins
            if (VisitsTotal(i).ne.0) then
               p_jnorm(i)=p_jnonmin(i)/(sum(p_jnonmin)*HistInt)
            else
               p_jnorm(i)=0.d0
            endif
        enddo


      open(unit=421, file='BL.Pjnorm.lnGj.Djnm.Djm.VT.his', status='unknown')
      do i=1, Hbins
         write(421, '(6G20.10)')  BinLabel(i), p_jnorm(i), lng_j(i), d_jnonmin(i), d_jmin(i), VisitsTotal(i)
      enddo
         close(421)

      end subroutine Record_Stat

!---======================================---
    subroutine SaveBinStructures(CurrentPointEnergy, CurrentPointCoordinates, BinIndex, MinimaNumber, WriteStruct)
        use commons, only:Natoms, hbins,BQMAX, CHRMMT, ZSYM
        implicit none

        integer binindex, jb,  MinimaNumber(Hbins), jc
        integer, save :: binensaved=0
        real(8) CurrentPointEnergy, CurrentPointCoordinates(3*Natoms, 1), Q4orderparam
        logical newbinenergy, yesno, WriteStruct
        character (len =256)  binname
        character (len= 10)  istr
        real(8),allocatable, save :: binenergies(:)
        integer,allocatable, save :: binenImportanceIndex(:)
        
         if (binensaved.eq.0) then
           allocate(binenergies(100000000))
           binenergies=0.0
           allocate(binenImportanceIndex(100000000))
           binenImportanceIndex=0
         endif
      
         if (binensaved.eq.size(binenergies)) then
            return
         endif

         inquire(file='binenergies', exist=yesno)

         newbinenergy=.true.
         do jb=1, binensaved
             if (abs(binenergies(jb)-CurrentPointEnergy).le.BQMAX) then
                newbinenergy=.false.
             endif
         enddo

         if (newbinenergy) then
            binensaved=binensaved+1
            !call ORDERQ4(Natoms, CurrentPointCoordinates, Q4orderparam)
            binenergies(binensaved)=CurrentPointEnergy
            if (WriteStruct)  then
                if ((binensaved.eq.1).and.(yesno)) then
                   open(unit=1979,file='binenergies',status='old')
                else
                   open(unit=1979,file='binenergies',status='unknown',position='append')
                endif
                write(1979,'(2G20.10)') CurrentPointEnergy, Q4orderparam
                close (1979)
                write (istr, '(i10)') binindex
                binname="binstructures."//trim(adjustl(istr))
                if ((binensaved.eq.1).and.(yesno)) then
                   open(unit=1979,file=binname, status="old", form="formatted")
                else
                   open(unit=1979,file=binname, status="unknown", form="formatted", position="append")
                endif
                write(1979,'(I10)') natoms
                write(1979,'(A,1G20.10)') '     Structure with energy ', CurrentPointEnergy
                IF (CHRMMT) THEN
                   DO JB=1,NATOMS
                      write(1979,'(A,1X,3F20.10)') ZSYM(JB)(1:1),(CurrentPointCoordinates(3*(JB-1)+JC,1),JC=1,3)
                   ENDDO
                ELSE
                   write(1979,30) (CurrentPointCoordinates(JB,1),JB=1,3*NATOMS)
30                 format('LA ',3F20.10)
                ENDIF
                close (1979)
             endif
             MinimaNumber(BinIndex)=MinimaNumber(BinIndex)+1
         endif

! Importance Index 

 	 if ((.not.newbinenergy).and.(WriteStruct)) then 
		open(unit=1979,file="BinEnImportanceIndex", status="unknown", form="formatted")
		do jb=1, binensaved
	             if (abs(binenergies(jb)-CurrentPointEnergy).le.BQMAX) then
		     binenImportanceIndex(jb)=binenImportanceIndex(jb)+1
             	     endif
         	enddo
		do jb=1, binensaved
		write(1979,'(1G20.10, I10)') binenergies(jb), binenImportanceIndex(jb)
		enddo
		close (1979)
	 endif
 
         if (binensaved.eq.size(binenergies)) then
            deallocate(binenergies)
         endif

         end subroutine savebinstructures



      end subroutine BasinSampling

