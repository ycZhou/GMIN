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
      SUBROUTINE LJPARAMMF
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592654)
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      open(unit=15,file='mf.dat',status='old')

      READ(15,*)H
      READ(15,*)M
      READ(15,*)GD
      READ(15,*)L
      READ(15,*)A
      READ(15,*)B
      READ(15,*)T
      READ(15,*)EPSILON
      READ(15,*)RHOHAT
      READ(15,*)LAMBDA
      READ(15,*)MU

      FACT=MU*PI*RHOHAT

      CLOSE(15)

      RETURN
      END
C
C  Energy and Gradient for the Mean Field Potential.
C  
C
      SUBROUTINE MF(X,V,EMF,GTEST)
      USE commons
      IMPLICIT NONE

      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      LOGICAL GTEST
      INTEGER J1, J3, IJ
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS)
      DOUBLE PRECISION EMF,RMIN,RMAX
      DOUBLE PRECISION R(NATOMS),RR,RC,SUM,DD
      DOUBLE PRECISION VMF,D_R_VMF,D_RC_VMF
      DOUBLE PRECISION INT1,INT2,INT3,INT4,INT5
      DOUBLE PRECISION INFTY
      PARAMETER (INFTY=1.0D30)
      EXTERNAL FUNC1,FUNC2,FUNC3,FUNC4

      COMMON /RRR/RR,RC
      COMMON /IIJ/IJ


        RMAX=0.0D0
        EMF =0.0D0

         DO J1=1,NATOMS

           J3=3*J1
           DIST =DSQRT(X(J3-2)**2+X(J3-1)**2+X(J3)**2)
           R(J1)=DIST

           IF(DIST.GT.RMAX) THEN
             RMAX=DIST
             IJ  =J1
           ENDIF

         ENDDO

          RC    = RMAX+LAMBDA
          RADIUS= RC*RC
          SUM   = 0.0D0

         DO J1=1,NATOMS
           
              RR  =R(J1)

              RMIN=RC-RR
              RMAX=RC+RR

C CALCULATION OF POTENTIAL AND DERIVATIVES

             call qromb(func1,RMIN,RMAX,INT1)
             call qromb(func2,RMIN,RMAX,INT2)
             call qromb(func3,RMIN,RMAX,INT4)
             call qromb(func4,RMIN,RMAX,INT5)

C CALCULATION OF THE INTEGRAL I3 (TO INFINITY)

             call qromo(func1,RMAX,INFTY,INT3)

             VMF      = FACT * (INT1-INT2+2.0D0*INT3)
             
             D_R_VMF  = FACT * INT4
             D_RC_VMF = -FACT * RC * INT5 / RR 

             EMF   = EMF+VMF

             VT(J1)= VT(J1)+VMF/4.0D0

             DD  = D_R_VMF/RR

             J3     =3*J1 
             V(J3)  =DD*X(J3)  
             V(J3-1)=DD*X(J3-1)
             V(J3-2)=DD*X(J3-2)

             SUM = SUM + D_RC_VMF
             
          ENDDO
          
C************************
C PARTICULAR CASE RR=RMAX
C************************

      J3=3*IJ
      RR=R(IJ)
      DD=SUM/RR
      V(J3)  =V(J3)  +DD*X(J3)  
      V(J3-1)=V(J3-1)+DD*X(J3-1)
      V(J3-2)=V(J3-2)+DD*X(J3-2)

      call LJ_MF_DIFF(X)

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine LJDIFF calculates the cartesian second
C  derivative matrix analytically. Reduced units.
C
C*************************************************************************
C
      SUBROUTINE LJ_MF_DIFF(X)
      USE commons
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, IJ
      DOUBLE PRECISION X(3*NATOMS),R0(NATOMS), 
     1                 R2(NATOMS,NATOMS), 
     2                 R8(NATOMS,NATOMS), G(NATOMS,NATOMS),
     3                 R14(NATOMS,NATOMS), F(NATOMS,NATOMS)
      DOUBLE PRECISION INT4,INT5,INT6
      DOUBLE PRECISION D_R_VMF,D_RC_VMF,D_RC_RC_VMF,D_R_R_VMF,D_R_RC_VMF
      DOUBLE PRECISION D_R(NATOMS),D_RC(NATOMS)
      DOUBLE PRECISION D_R_R(NATOMS),D_RC_RC(NATOMS),D_R_RC(NATOMS)
      DOUBLE PRECISION SUM1,SUM2
      COMMON /D_MF/SUM1,SUM2
      DOUBLE PRECISION RR,RC,RMIN,RMAX
      COMMON /RRR/RR,RC
      EXTERNAL FUNC3,FUNC4,FUNC6
      double precision FUNC5
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

      N=NATOMS
C 
C  Store distance matrices.
C
      SUM1=0.D0
      SUM2=0.D0

      DO J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         J3=3*J1
         R0(J1)=(X(J3-2)**2+X(J3-1)**2+X(J3)**2)
C
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R8(J2,J1)=R2(J2,J1)**4
            R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
            R2(J1,J2)=R2(J2,J1)
         ENDDO
C
C  Calculate MF derivatives
C           
              RR  =SQRT(R0(J1))

              RMIN=RC-RR
              RMAX=RC+RR

             call qromb(func3,RMIN,RMAX,INT4)
             call qromb(func4,RMIN,RMAX,INT5)
             call qromb(func6,RMIN,RMAX,INT6)

             D_R_VMF  = FACT * INT4
             D_RC_VMF = -FACT * RC * INT5 / RR 
             D_R_RC_VMF = -FACT*RC/RR*
     1       (RMAX*FUNC5(RMAX)+RMIN*FUNC5(RMIN)-INT5/RR)
             D_RC_RC_VMF =-FACT/RR*
     1       (RC*RMAX*FUNC5(RMAX)-RC*RMIN*FUNC5(RMIN)+INT5)
             D_R_R_VMF =FACT/RR*
     1       (-RC*RMAX*FUNC5(RMAX)+RC*RMIN*FUNC5(RMIN)+INT6/RR**2)
C
             D_R    (J1)=D_R_VMF
             D_RC   (J1)=D_RC_VMF
             D_RC_RC(J1)=D_RC_RC_VMF
             D_R_RC (J1)=D_R_RC_VMF
             D_R_R  (J1)=D_R_R_VMF

             SUM1=SUM1+D_RC_RC_VMF
             SUM2=SUM2+D_RC_VMF

      ENDDO

C     CALL LJ_MF_S(G,F,R0,R2,R14,R8,X)

      RETURN
      END

C*****************************************************************************

C     SUBROUTINE LJ_MF_S(G,F,R0,R2,R14,R8,X)
C     USE commons
C     IMPLICIT NONE
C     INTEGER N, J1, J2, J3, J4, J5, J6, IJ, II, JJ
C     DOUBLE PRECISION G(NATOMS,NATOMS), R14(NATOMS,NATOMS), R8(NATOMS,NATOMS),
C    1                 R2(NATOMS,NATOMS), F(NATOMS,NATOMS), R0(NATOMS),
C    2                 X(3*NATOMS),DUMMYLJ,DUMMYMF
C     DOUBLE PRECISION D_R(NATOMS),D_RC(NATOMS)
C     DOUBLE PRECISION D_R_R(NATOMS),D_RC_RC(NATOMS),D_R_RC(NATOMS)
C     DOUBLE PRECISION SUM1,SUM2
C     COMMON /D_MF/D_R,D_RC,D_R_R,D_RC_RC,D_R_RC,SUM1,SUM2
C     DOUBLE PRECISION DIST,EELJ,EEMF,DUMMY,EPP,EPM,EMP,EMM,DX,R6,DDNUM
C     DOUBLE PRECISION EMF,RMIN,RMAX,DUMMYMF1,DUMMYMF2
C     DOUBLE PRECISION R(NATOMS),RR,RC,SUM,DD
C     DOUBLE PRECISION VMF,D_R_VMF,D_RC_VMF
C     DOUBLE PRECISION INT1,INT2,INT3,INT4,INT5,INT6
C     DOUBLE PRECISION INFTY
C     COMMON /RRR/RR,RC
C     COMMON /IIJ/IJ
C     PARAMETER (INFTY=1.0D30)
C     EXTERNAL FUNC1,FUNC2,FUNC3,FUNC4
C     DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
C     COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU
C
C     N=NATOMS
C     DO J1=1,N
C        G(J1,J1)=0.0D0
C        F(J1,J1)=0.0D0
C        DO J2=J1+1,N 
C           F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
C           F(J1,J2)=F(J2,J1)
C           G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
C           G(J1,J2)=G(J2,J1)
C        ENDDO
C     ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
C    

C     DO J1=1,N
C        DO J2=1,3
C    J3=3*(J1-1)+J2
C    DUMMYLJ=0.0D0          
C    DO J4=1,N
C              DUMMYLJ=DUMMYLJ+F(J4,J1)*R2(J4,J1)*
C    1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
C           ENDDO

C******************************************
C     DUMMYMF=0.0D0
C     IF(J1.EQ.IJ) THEN
C      DUMMYMF=(D_R_R(J1)+2.0D0*D_R_RC(J1)+SUM1)*X(J3)*X(J3)/R0(J1)
C    2  +(D_R(J1)+SUM2)*(1.0D0-X(J3)*X(J3)/R0(J1))/SQRT(R0(J1))
C     ELSE
C     DUMMYMF=D_R_R(J1)*X(J3)*X(J3)/R0(J1)
C    2  +D_R(J1)*(1.0D0-X(J3)*X(J3)/R0(J1))/SQRT(R0(J1))
C     ENDIF
C******************************************
C         HESS(J3,J3)=DUMMYLJ+DUMMYMF
C       ENDDO
C     ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
C     DO J1=1,N
C DO J2=1,3
C    J3=3*(J1-1)+J2
C    DO J4=J2+1,3
C       DUMMYLJ=0.0D0
C       DO J5=1,N
C                 DUMMYLJ=DUMMYLJ + F(J5,J1)*R2(J5,J1)* 
C    1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
C              ENDDO
C              
C******************************************
C     DUMMYMF=0.0D0
C     II=3*(J1-1)+J4
C     IF(J1.EQ.IJ) THEN
C      DUMMYMF=(D_R_R(J1)+2.0D0*D_R_RC(J1)+SUM1)*X(II)*X(J3)/R0(J1)
C    2  -(D_R(J1)+SUM2)*X(II)*X(J3)/R0(J1)/SQRT(R0(J1))
C     ELSE
C     DUMMYMF=D_R_R(J1)*X(II)*X(J3)/R0(J1)
C    2  -D_R(J1)*X(II)*X(J3)/R0(J1)/SQRT(R0(J1))
C     ENDIF
C******************************************
C     HESS(3*(J1-1)+J4,J3)=DUMMYLJ+DUMMYMF
C           ENDDO
C        ENDDO
C     ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C

C     DO J1=1,N
C DO J2=1,3
C    J3=3*(J1-1)+J2
C    DO J4=J1+1,N
C              II=3*(J4-1)+J2
C              HESS(II,J3)=-F(J4,J1)*R2(J4,J1)*
C    1                           (X(J3)-X(II))**2-G(J4,J1) 

C******************************************
C     DUMMYMF=0.0D0
C     IF(J1.EQ.IJ.OR.J4.EQ.IJ) THEN
C     DUMMYMF=D_R_RC(J1)*X(II)*X(J3)/SQRT(R0(J1))/SQRT(R0(J4))
C     ELSE
C     DUMMYMF=0.0D0
C     ENDIF
C      HESS(II,J3)=HESS(II,J3)+DUMMYMF
C******************************************
C           ENDDO
C        ENDDO
C     ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
C     DO J1=1,N
C        DO J2=1,3
C           J3=3*(J1-1)+J2
C           DO J4=J1+1,N
C              DO J5=1,J2-1
C                 J6=3*(J4-1)+J5
C                 HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
C    1                    *(X(J3)-X(3*(J4-1)+J2))
C    2                    *(X(3*(J1-1)+J5)-X(J6))
C                 HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)

C******************************************
C     DUMMYMF1=0.0D0
C     DUMMYMF2=0.0D0
C     IF(J1.EQ.IJ.OR.J4.EQ.IJ) THEN
C     DUMMYMF1=D_R_RC(J1)*X(J6)*X(J3)/SQRT(R0(J1))/SQRT(R0(J4))
C     DUMMYMF2=D_R_RC(J1)*X(3*(J4-1)+J2)*X(3*(J1-1)+J5)/SQRT(R0(J1))/SQRT(R0(J4))
C     ELSE
C     DUMMYMF1=0.0D0
C     DUMMYMF2=0.0D0
C     ENDIF

C     HESS(J6,J3)=HESS(J6,J3)+DUMMYMF1
C     HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(3*(J4-1)+J2,3*(J1-1)+J5)+DUMMYMF2
C******************************************
C              ENDDO
C           ENDDO
C        ENDDO
C     ENDDO
C
C  Symmetrise Hessian
C
C     DO J1=1,3*N
C        DO J2=J1+1,3*N
C           HESS(J1,J2)=HESS(J2,J1)
C        ENDDO
C     ENDDO

C     RETURN
C     END

c***************************************
      DOUBLE PRECISION FUNCTION func1(y)
      USE commons
      IMPLICIT NONE
c***************************************
      DOUBLE PRECISION Y,Y2,GDR,VLJ
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU
      
        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

       func1   =Y2*GDR*VLJ
                     
      END

c***************************************
      DOUBLE PRECISION FUNCTION func2(y)
      USE commons
      IMPLICIT NONE
c***************************************

      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      INTEGER IJ
      COMMON /RRR/RR,RC
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ

        func2 = Y2*FF*(RC*RC-RR*RR-Y2)/(2.0D0*RR*Y)

        END

c***************************************
      DOUBLE PRECISION FUNCTION func3(y)
      USE commons
      IMPLICIT NONE
c***************************************

      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      INTEGER IJ
      COMMON /RRR/RR,RC
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ

        func3 = Y2*FF*(RC*RC+RR*RR-Y2)/(2.0D0*RR*RR*Y)

      END

c***************************************
      DOUBLE PRECISION FUNCTION func4(y)
      USE commons
      IMPLICIT NONE
c***************************************

      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF    =GDR*VLJ
        func4 = Y*FF

      END

c***************************************
      DOUBLE PRECISION FUNCTION func5(y)
      USE commons
      IMPLICIT NONE
c***************************************

      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF    =GDR*VLJ
        func5 =FF

      END

c***************************************
      DOUBLE PRECISION FUNCTION func6(y)
      USE commons
      IMPLICIT NONE
c***************************************

      DOUBLE PRECISION Y,Y2,GDR,VLJ,FF,RC,RR
      DOUBLE PRECISION t1,t2,t3,t12,t15,t10,t6,t4,cgret      
      INTEGER IJ
      COMMON /RRR/RR,RC
      DOUBLE PRECISION H,M,GD,L,A,T,B,EPSILON,LAMBDA,RHOHAT,FACT,MU
      COMMON /COEFF/ H,M,GD,L,A,T,EPSILON,RHOHAT,LAMBDA,B,FACT,MU

        Y2=Y*Y
        IF(Y.GT.H) THEN

        t2 = y / h        
        t3 = t2**m
        t10= t2 - 0.1D1
        t12= exp(-a * t10)
        t15= cos(b * t10)
        cgret = 0.1D1 + 0.10D1 / t3 * dble(gd - 1 - l)+
     1  (t2 - 0.1D1 + dble(l)) / y * h * t12 * t15
        GDR=cgret
        GOTO 11

      ENDIF
        
        t4 = (y / h - 0.1D1)**2        
        t6 = exp(-t * t4)
        cgret = gd * t6
        GDR = cgret

 11     CONTINUE
        
        t1 = y**2        
        t2 = t1**2
        t3 = t2**2
        cgret = 0.4D1 * epsilon * (0.10D1 / t3 / t2 - 0.10D1 / t2 / t1)
        VLJ=cgret

        FF   =GDR*VLJ

        func6 = Y*FF*(Y2-RC*RC)

        END

c************************************
      SUBROUTINE qromb(func,aa,bb,ss)
      IMPLICIT NONE
c************************************

      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION aa,bb,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.0D-9, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
CU    USES polint,trapzd
      INTEGER j
      DOUBLE PRECISION dss,hh(JMAXP),s(JMAXP)
      hh(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,aa,bb,s(j),j)
        if (j.ge.K) then
          call polint(hh(j-KM),s(j-KM),K,ZERO,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        hh(j+1)=0.25*hh(j)
11    continue
      PRINT*, 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.

c**************************************
      SUBROUTINE trapzd(func,aa,bb,s,n)
      IMPLICIT NONE
c**************************************

      INTEGER n
      DOUBLE PRECISION aa,bb,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,xx
      if (n.eq.1) then
        s=0.5*(bb-aa)*(func(aa)+func(bb))
      else
        it=2**(n-2)
        tnm=it
        del=(bb-aa)/tnm
        xx=aa+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(xx)
          xx=xx+del
11      continue
        s=0.5*(s+(bb-aa)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.

c****************************************
      SUBROUTINE polint(xa,ya,n,xx,yy,dy)
      IMPLICIT NONE
c****************************************

      INTEGER n,NMAX
      DOUBLE PRECISION dy,xx,yy,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(xx-xa(1))
      do 11 i=1,n
        dift=abs(xx-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      yy=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-xx
          hp=xa(i+m)-xx
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)PRINT*, 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        yy=yy+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.

c****************************************
      SUBROUTINE midinf(funk,inf,sup,s,n)
      IMPLICIT NONE
c****************************************

      INTEGER n
      DOUBLE PRECISION inf,sup,s,funk
      EXTERNAL funk
      INTEGER it,j
      DOUBLE PRECISION aa,bb,ddel,del,sum,tnm,func,xx
      func(xx)=funk(1./xx)/xx**2
      bb=1./inf
      aa=1./sup
      if (n.eq.1) then
        s=(bb-aa)*func(0.5*(aa+bb))
      else
        it=3**(n-2)
        tnm=it
        del=(bb-aa)/(3.*tnm)
        ddel=del+del
        xx=aa+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(xx)
          xx=xx+ddel
          sum=sum+func(xx)
          xx=xx+del
11      continue
        s=(s+(bb-aa)*sum/tnm)/3.
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.

c*******************************************
      SUBROUTINE qromo(func,aa,bb,ss)
      IMPLICIT NONE
c*******************************************

      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION aa,bb,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.0D-9, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
CU    USES polint
      INTEGER j
      DOUBLE PRECISION dss,hh(JMAXP),s(JMAXP)
      hh(1)=1.
      do 11 j=1,JMAX
        call midinf(func,aa,bb,s(j),j)
        if (j.ge.K) then
          call polint(hh(j-KM),s(j-KM),K,ZERO,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        hh(j+1)=hh(j)/9.
11    continue
      PRINT*, 'too many steps in qromo'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.





