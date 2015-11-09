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

MODULE nrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0d0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
        complex(dpc), parameter :: CZ=(0.0_sp,0.0_sp), CI=(0.0_sp,1.0_sp)
	TYPE sprs2_sp
		INTEGER(I4B) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
	TYPE sprs2_dp
		INTEGER(I4B) :: n,len
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE nrtype

SUBROUTINE DIFFRACT ( X, V, ELJ, GTEST, SECT, NATOMS )
use nrtype
IMPLICIT NONE

include 'NP.H'
integer          :: natoms 
LOGICAL                    :: GTEST,SECT
DOUBLE PRECISION          :: X(3*NATOMS), V(3*NATOMS)
real*8                     :: elj 
logical                    :: ONLYONCE = .true.
real*8                     :: kEX
common /ENE/                  kEX
real*8, dimension (ncoso,ncosi)    :: intEX
real*8, dimension (ncoso)          :: cosoEX 
real*8, dimension (ncosi)          :: cosiEX 

common /EXP/ cosoEX, cosiEX, intEX

!
real*8, dimension(3*NATOMS)      :: xv

INTERFACE

!  s = \sum_{i} f_{i} (\vec)
   FUNCTION funcs(xv,NATOMS)
      REAL*8, DIMENSION(:), INTENT(IN)       :: xv
      REAL*8                                 :: funcs
      INTEGER NATOMS
   END FUNCTION funcs

!  grad_{i} s = \frac { \partial s }{ \partial x_{i} }
   FUNCTION dfuncs(xv,NATOMS)
      REAL*8, DIMENSION(:), INTENT(IN)      :: xv
      REAL*8, DIMENSION(size(xv))           :: dfuncs
      INTEGER NATOMS
   END FUNCTION dfuncs

END INTERFACE

!  Make pseudo-experimental database (only once).
if ( ONLYONCE ) then
      call makeEX(NATOMS)
      ONLYONCE = .false.
endif

xv  =         x
elj =  funcs (xv, NATOMS)
v   = dfuncs (xv, NATOMS)

return

end subroutine DIFFRACT

FUNCTION funcs ( xv, NATOMS )
! The Cost Function
! f ( \vec x )
IMPLICIT NONE
include 'NP.H'

integer                        :: io, ii, NATOMS
REAL*8, DIMENSION(:), INTENT(IN) :: xv
REAL*8                           :: funcs

complex*16                       :: fsc
complex*16,dimension(3*NATOMS)         :: dfscCJ
real*8, parameter                :: favg = 1.e6
REAL*8                           :: coso, cosi 

real*8                           :: kEX
common /ENE/ kEX
real*8, dimension (ncoso)        :: cosoEX
real*8, dimension (ncosi)        :: cosiEX
real*8, dimension (ncoso,ncosi)  :: intEX
common /EXP/ cosoEX, cosiEX, intEX

funcs = 0.0
do ii = 1, ncosi
   cosi = cosiEX (ii) 
   do io = 1, ncoso
      coso = cosoEX (io) 
      call msmott_angle ( 3*NATOMS, 3*NATOMS-3, kEX, cosi, coso, xv, fsc )
      funcs = funcs + ( ( Abs (fsc) )**2  - intEX(io,ii)  )**2
   enddo
enddo

END FUNCTION funcs
FUNCTION dfuncs ( xv, NATOMS )
use nrtype
IMPLICIT NONE
include 'NP.H'
INTEGER NATOMS

INTERFACE
!  s = \sum_{i} f_{i} (\vec)
   FUNCTION funcs(xv,NATOMS)
      REAL*8, DIMENSION(:), INTENT(IN)       :: xv
      REAL*8                                 :: funcs
   END FUNCTION funcs
END INTERFACE

REAL*8, DIMENSION(:), INTENT(IN)    :: xv
REAL*8, DIMENSION(size(xv))         :: dfuncs

complex*16                          :: fsc
complex*16 , dimension (3*NATOMS)         :: dfsc
REAL*8                              :: coso, cosi
real*8, parameter                   :: favg = 1.e12
integer                          :: ix, io, ii

real*8                           :: kEX
common /ENE/ kEX
real*8, dimension (ncoso)        :: cosoEX
real*8, dimension (ncosi)        :: cosiEX
real*8, dimension (ncoso,ncosi)  :: intEX
common /EXP/ cosoEX, cosiEX, intEX

do ix = 1, size(xv)

   dfuncs (ix) = 0.0
   do ii = 1, ncosi
      cosi = cosiEX (ii) 
      do io = 1, ncoso
         coso = cosoEX (io) 

         call msmott_angle (  3*NATOMS, 3*NATOMS-3,  kEX, cosi, coso, xv, fsc)
         call dmsmott_angle ( 3*NATOMS, 3*NATOMS-3,  kEX, cosi, coso, xv, dfsc)

         dfuncs (ix) =  dfuncs (ix) + &
                        2. * ( ( Abs (fsc) )**2  - intEX(io,ii)  ) * &
                   Real( fsc*CONJG(dfsc(ix)) + CONJG(fsc)*dfsc(ix))

      enddo   ! io
   enddo   ! ii

enddo  ! ix


! Parabolas para los desfasajes (1:3*NATOMS-2), 
! distancia (3*NATOMS-1) y angulo (3*NATOMS), en el intervalo (0,Pi).

!do ix = 1, 3*NATOMS
!   if ( xv(ix) .gt.  PI ) then
!      dfuncs (ix) = dfuncs (ix) * ( 1. + favg * ( xv(ix) - PI )**2 ) + &
!                    2.*favg* ( xv(ix) - PI ) * funcs (xv)
!   endif
!   if ( xv(ix) .lt.  0. ) then
!      dfuncs (ix) = dfuncs (ix) * ( 1. + favg * ( xv(ix) + 0. )**2 ) + &
!                    2.*favg* ( xv(ix) + 0. ) * funcs (xv)
!   endif
!enddo

END FUNCTION dfuncs
SUBROUTINE msmott_angle  ( NP, lmax, kEX, cosi, coso, xvec, amplitude)

USE nrtype
! USE ncalls
IMPLICIT NONE

! Double precision for David Wales.
INTEGER , INTENT(IN)    :: NP, lmax
REAL*8 , INTENT(IN) , DIMENSION( NP )    :: xvec
REAL*8 , INTENT(IN)        :: cosi, coso, kEX
COMPLEX*16 , INTENT(OUT) :: amplitude

! Internal
INTEGER                      :: il
REAL(SP)                     :: kr, ith, r0, energia, Voi, alpha, cosa
REAL(SP)                     :: cosoi, cosoa, cosia, sinia, sinoa, sina
REAL(SP), DIMENSION (0:lmax) :: LegendrePo, LegendrePi, LegendrePoi, delta
COMPLEX(DPC)                 :: til, amplth, eikrr0, eikrcoso, &
                                eikrcosi, eikrcosoi, kExc
COMPLEX(DPC)                 :: fo, fi, f_o, f_i, f_1, foi, dumi, &
                                dtil1, dtilo, dtili, corchete, den, &
                                parentesis
COMPLEX(DPC)                 :: C1=(1.,0.)

!------------------------------------------------------------

r0 = xvec(NP-1)
alpha = xvec(NP) - 0.5*pi           !!!!!!!!!
cosa = COS(alpha)
sina = SIN(alpha)

! cosoi = cos(theta_in-theta_out)
cosoi = coso*cosi+SQRT( (1.-coso*coso ) * (1.-cosi*cosi ) )
kr = kEX * r0
!eikrr0 = EXP ( CI * kr ) / r0
!   Introduzco Voi
    energia = kEX * kEX * 0.5
    Voi = 0.0
    kExc = SQRT( CMPLX( 2.0*energia, -2.0*Voi))
!print *, 'r0', r0
    eikrr0 = EXP ( CI * kExc * r0 ) / r0

! cosoa = cos(theta_out + alpha)
! cosia = cos(theta_in + alpha)
cosoa = coso*cosa - SQRT( 1.-coso*coso ) * sina 
cosia = cosi*cosa - SQRT( 1.-cosi*cosi ) * sina 
eikrcoso = EXP ( CI * kr * cosoa)
eikrcosi = EXP ( CI * kr * cosia)
eikrcosoi = CONJG(eikrcoso) * eikrcosi

! sinoa = sin(theta_out + alpha)
! sinia = sin(theta_in + alpha)
sinoa = SQRT( 1.-coso*coso )*cosa + coso * sina
sinia = SQRT( 1.-cosi*cosi )*cosa + cosi * sina

! Calculo los factores de scattering atomicos:
! f( cosoa ) == fo, f (cosia) == fi, f(-cosoa) == f_o, f(-cosia) == f_i, 
! f(cos(theta_in+theta_out)) == foi, f(R,-R) == f_1

CALL Plegendre ( cosoa, LegendrePo, lmax )
CALL Plegendre ( cosia, LegendrePi, lmax )
CALL Plegendre ( cosoi, LegendrePoi, lmax )

fo = CZ
fi = CZ
f_o = CZ
f_i = CZ
foi = CZ
f_1 = CZ
delta(0:lmax) = xvec(1:NP-1)
DO il = 0, lmax
   til = -CI * 0.5 * (EXP( CI * 2. * delta (il) ) - C1) * ( 2.*il + 1.) / kEX
   fo = fo + til * LegendrePo(il)
   f_o = f_o + til * LegendrePo(il) * (-1)**il
   fi = fi + til * LegendrePi(il)
   f_i = f_i + til * LegendrePi(il) * (-1)**il
   foi = foi + til * LegendrePoi(il)
   f_1 = f_1 + til * (-1)**il
END DO

! Parentesis
parentesis = eikrcosoi*fo*fi*f_1+f_o*fi*f_1

! Corchete
corchete =  CONJG(eikrcoso)*fo*fi + eikrcosi*f_o*f_i + &
       eikrr0 * parentesis

! Denominador "1-X"
!WRITE(*,'(a,e14.4)') 'r0',r0
!WRITE(*,'(a,2e14.4)') 'eikrr0',eikrr0
!WRITE(*,'(a,2e14.4)') 'f_1',f_1
den = C1 - (eikrr0*f_1)**2

! Total scattered amplitude
amplitude = foi * ( C1 + eikrcosoi ) + eikrr0*corchete / den

RETURN
END SUBROUTINE msmott_angle
SUBROUTINE dmsmott_angle  ( NP, lmax, kEX, cosi, coso, xvec, dvec)
USE nrtype
IMPLICIT NONE

! Double precision for David Wales.
INTEGER , INTENT(IN)    :: NP, lmax
REAL*8 , INTENT(IN) , DIMENSION( NP )    :: xvec
REAL*8 , INTENT(IN)        :: cosi, coso, kEX
COMPLEX*16 , DIMENSION( NP ) ,INTENT(OUT) :: dvec

! Internal
INTEGER                      :: il
REAL(SP)                     :: kr, ith, r0, energia, Voi, cosa, alpha
REAL(SP)                     :: cosoi, cosoa, cosia, sinia, sinoa, sina
REAL(SP), DIMENSION (0:lmax) :: LegendrePo, LegendrePi, LegendrePoi, delta, &
                                dLegendrePo, dLegendrePi
COMPLEX(DPC)                 :: til, amplth, eikrr0, eikrcoso, eikrcosi, &
                                eikrcosoi, kEXc
COMPLEX(DPC)                 :: fo, fi, f_o, f_i, f_1, foi, dumi, &
                                dtil1, dtilo, dtili, corchete, den, parentesis
COMPLEX(DPC)                 :: dfo, dfi, df_o, df_i 
COMPLEX(DPC)                 :: C1=(1.,0.)

!------------------------------------------------------------

r0 = xvec(NP-1)
alpha = xvec(NP) - pi * 0.5               !!!!!!!!!!
cosa = COS(alpha)
sina = SIN(alpha)
! cosoi = cos(theta_in-theta_out)
cosoi = coso*cosi+SQRT( (1.-coso*coso ) * (1.-cosi*cosi ) )
kr = kEX * r0
!eikrr0 = EXP ( CI * kr ) / r0
!   Introduzco Voi
    energia = kEX * kEX * 0.5
    Voi = 0.0
    kExc = SQRT( CMPLX( 2.0*energia, -2.0*Voi))
    eikrr0 = EXP ( CI * kExc * r0 ) / r0

! cosoa = cos(theta_out + alpha)
! cosia = cos(theta_in + alpha)
cosoa = coso*cosa - SQRT(1.-coso*coso) * sina
cosia = cosi*cosa - SQRT(1.-cosi*cosi) * sina
eikrcoso = EXP ( CI * kr * cosoa)
eikrcosi = EXP ( CI * kr * cosia)
eikrcosoi = CONJG(eikrcoso) * eikrcosi

! sinoa = sin(theta_out + alpha)
! sinia = sin(theta_in + alpha)
sinoa = SQRT( 1.-coso*coso )*cosa + coso * sina
sinia = SQRT( 1.-cosi*cosi )*cosa + cosi * sina

! Calculo los factores de scattering atomicos:
! f( cosoa ) == fo, f (cosia) == fi, f(-cosoa) == f_o, f(-cosia) == f_i, 
! f(cos(theta_in+theta_out)) == foi, f(R,-R) == f_1

CALL Plegendre ( cosoa, LegendrePo, lmax )
CALL Plegendre ( cosia, LegendrePi, lmax )
CALL Plegendre ( cosoi, LegendrePoi, lmax )

fo = CZ
fi = CZ
f_o = CZ
f_i = CZ
foi = CZ
f_1 = CZ
delta(0:lmax) = xvec(1:NP-1)
DO il = 0, lmax
   til = -CI * 0.5 * (EXP( CI * 2. * delta (il) ) - C1) * ( 2.*il + 1.) / kEX
   fo = fo + til * LegendrePo(il)
   f_o = f_o + til * LegendrePo(il) * (-1)**il
   fi = fi + til * LegendrePi(il)
   f_i = f_i + til * LegendrePi(il) * (-1)**il
   foi = foi + til * LegendrePoi(il)
   f_1 = f_1 + til * (-1)**il
END DO

! Parentesis
parentesis = eikrcosoi*fo*fi*f_1+f_o*fi*f_1

! Corchete
corchete =  CONJG(eikrcoso)*fo*fi + eikrcosi*f_o*f_i + &
       eikrr0*parentesis

! Denominador "1-X"
den = C1 - (eikrr0*f_1)**2

! Loop over deltas to generate the gradient
DO il = 0, lmax

   dtil1 = EXP( CI * 2. * delta (il) ) * ( 2.*il + 1.) / kEX
   dtilo = dtil1 * LegendrePo(il) 
   dtili = dtil1 * LegendrePi(il) 

   dvec(il+1) = dtil1 * LegendrePoi(il) * ( C1 + eikrcosoi)

   dumi = CONJG(eikrcoso)*( dtili*fo + fi*dtilo ) + &
          eikrcosi * ( (-1)**il*dtili*f_o + f_i*dtilo*(-1)**il ) + &
          eikrr0 * ( eikrcosoi * ( dtilo*fi*f_1 + fo*dtili*f_1 + fo*fi*dtil1*(-1)**il ) &
                     + dtilo*fi*f_1*(-1)**il + f_o*dtili*f_1 + f_o*fi*dtil1*(-1)**il  )
   dumi = dumi * eikrr0 / den
   dvec(il+1) = dvec(il+1) + dumi

   dumi = 2.0 * f_1 * dtil1*(-1)**il * eikrr0**3 * corchete / den**2
   dvec(il+1) = dvec(il+1) + dumi

END DO

! Derivada respecto de R0
dvec(lmax+2) = eikrcosoi*foi*CI*kEX*(cosia-cosoa)

dumi = -CI*kEX*cosoa*CONJG(eikrcoso)*fo*fi + CI*kEX*cosia*eikrcosi*f_o*f_i + &
       eikrr0*CI*kEX*(cosia-cosoa)*eikrcosoi*fo*fi*f_1 + &
       eikrr0 * (CI*kEXc*r0-C1)*parentesis / r0 
dumi = dumi * eikrr0 / den
dvec(lmax+2) = dvec(lmax+2) + dumi

dumi = eikrr0 * (CI*kEXc*r0-C1)*corchete / ( r0 * den )
dvec(lmax+2) = dvec(lmax+2) + dumi

dumi = 2.0 * (f_1/den)**2 * eikrr0**3 * (CI*kEXc*r0-C1)*corchete / r0
dvec(lmax+2) = dvec(lmax+2) + dumi

! Derivada respecto de alpha ( dvec(lmax+3) )
! Primero necesito las derivadas de los factores de scattering
CALL dPlegendre ( cosoa, dLegendrePo, lmax)
CALL dPlegendre ( cosia, dLegendrePi, lmax)
dfo = CZ
dfi = CZ
df_o = CZ
df_i = CZ
DO il = 0, lmax
   delta (il) = xvec (il+1)
   til = -CI * 0.5 * (EXP( CI * 2. * delta (il) ) - C1) * ( 2.*il + 1.) / kEX
   dfo = dfo - sinoa * til * dLegendrePo(il)
   df_o = df_o - sinoa * til * dLegendrePo(il) * (-1)**il
   dfi = dfi - sinia * til * dLegendrePi(il)
   df_i = df_i - sinia * til * dLegendrePi(il) * (-1)**il
END DO

dvec(lmax+3) = foi * CI * kr * ( -sinia+sinoa) * eikrcosoi 

dumi = CONJG(eikrcoso)*( CI*kr*sinoa*fi*fo + dfi*fo + fi*dfo ) + &
       eikrcosi*( -CI*kr*sinia*f_i*f_o + df_i*f_o + f_i*df_o ) + &
       eikrr0 * f_1 * (eikrcosoi* (dfi*fo+fi*dfo+CI*kr*(-sinia+sinoa)*fi*fo) + &
                       dfi*f_o+fi*df_o ) 
dumi = dumi * eikrr0 / den
dvec(lmax+3) = dvec(lmax+3) + dumi

RETURN
END SUBROUTINE dmsmott_angle

subroutine makeEX(NATOMS)
IMPLICIT NONE

!        Especifica el numero de parametros: NP, 
!        y el numero de datos: ncoso, ncosi
include 'NP.H'

! Valor minimo de global en presencia de ruido (0 si ruido es 0). 
real*8                               :: minimorum 
INTEGER NATOMS

real*8, dimension (3*NATOMS) :: xvEX
INTEGER                            :: i

! Datos experimentales, etc
complex*16                         :: fsc
complex*16,dimension(3*NATOMS)           :: dfscCJ
real*8                             :: PI=3.14159
real*8                             :: ruido = 0.0
real*8, dimension (ncoso,ncosi)    :: intEX
real*8, dimension (ncoso)          :: cosoEX 
real*8, dimension (ncosi)          :: cosiEX 
real*8                             :: coso, cosi 
integer                            :: ii, io
integer                            :: idum 
real*8                             :: random1, dummy
external                                random1

real*8                             :: kEX
!                         Shared with MAIN, funcs, dfuncs
common /ENE/ kEX
common /EXP/ cosoEX, cosiEX, intEX

   kEX = 1.
   open(unit=1,file='xvEX.inp',status='old')
      read(1,*) xvEX
   close(1)

      write(6,"(/'Making pseudo-experimental data base...')")
      write(6,"(/'Deltas-EX (xv)  ', 6(f8.4,1x)/)") xvEX

   ! Genera datos experimentales simulados
   do ii = 1, ncosi
      cosiEX (ii) = -1.
   enddo 

   dummy = -1.
   do io = 1, ncoso
      dummy = dummy + 1./float(ncoso)
      cosoEX (io) = dummy
   enddo

   minimorum = 0.0
   do ii = 1, ncosi
      cosi = cosiEX (ii)
      do io = 1, ncoso
         coso = cosoEX (io)
         call msmott_angle ( 3*NATOMS, 3*NATOMS-3, kEX, cosi, coso, xvEX, fsc )
!        call fsc1TH ( kEX, coso, xvEX, fsc, dfscCJ )
         intEX (io,ii) = ( Abs (fsc) )**2 
!                      * ( 1. + ruido * (random1(idum)-0.5) ) 
         minimorum = minimorum + ( intEX (io,ii) - ( Abs (fsc) )**2 )**2
      enddo
   enddo

return

end subroutine makeEX

subroutine PLegendre ( cosx, LegendreP, lmax )
! Calcula polinomios de Legendre de argumento real.
! La recurrencia evita errores de redondedo excesivos.
use nrtype
implicit none
integer(i4b)                   :: lmax, il
real(sp)                       :: cosx
real(sp), dimension (0:lmax)   :: LegendreP 
LegendreP (0) = 1.0
if ( lmax > 0 ) LegendreP (1) = cosx
if ( lmax > 1 ) then
  do il = 1, lmax-1
     LegendreP (il+1) =  2. * cosx * LegendreP (il) -   &
                                   LegendreP (il-1) - &
               ( cosx*LegendreP(il) - LegendreP(il-1) ) / float(il+1)
  enddo
end if
return
end
!--------------------------------------------------------
subroutine dPLegendre ( cosx, dLegendreP, lmax )
! Calcula derivadas de polinomios de Legendre de argumento real.
! La recurrencia evita errores de redondedo excesivos.
use nrtype
implicit none
integer(i4b)                   :: lmax, il
real(sp)                       :: cosx
real(sp), dimension (0:lmax)   :: dLegendreP, LegendreP

dLegendreP (0) = 0.0
IF ( lmax > 0 ) dLegendreP (1) = 1.0

IF ( lmax > 1 ) THEN
CALL PLegendre ( cosx, LegendreP, lmax )
do il = 2, lmax
   dLegendreP (il) = cosx * dLegendreP (il-1) + il * LegendreP(il-1)
enddo
END IF

return
end
