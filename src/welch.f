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
C
C  Energy and derivatives of the binary salt potential described
C  by Welch et al, JCP, 94, 4980, 1976 and Phillips et al, JCP,
C  94, 4980, 1991.
C  Energy is in hartree, length in Bohr.
C
      SUBROUTINE WEL(X, V, POTEL, GTEST, STEST)
      USE commons
      IMPLICIT NONE 
      INTEGER N, I, J
C             , J1, K, L, J2, J3
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), POTEL, XQ(NATOMS), ALPHA(NATOMS),
     1                 AC(NATOMS,NATOMS), Q(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS),
     2                 XMU(3*NATOMS), XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS),
     3                 RRD5(NATOMS,NATOMS)
C                      , DIF, TEMP1, V1, V2, XMU1(3*NATOMS), XMU2(3*NATOMS)
      LOGICAL GTEST, STEST

      N=NATOMS
      DO I=1,N
         IF (ZSYM(I).EQ.'MI') THEN
            Q(I)=-1.0D0
         ELSE IF (ZSYM(I).EQ.'PL') THEN
            Q(I)=1.0D0
         ELSE
            WRITE(*,'(A)') ' All atoms must be type PL or MI for Welch'
            STOP
         ENDIF
      ENDDO

      DO I=1,N
         IF (ZSYM(I).EQ.'PL') THEN 
            XQ(I)=XQP
            ALPHA(I)=ALPHAP
         ELSE
            XQ(I)=XQM
            ALPHA(I)=ALPHAM
         ENDIF
         DO J=I,N
            IF (ZSYM(I).EQ.'PL') THEN 
               IF (ZSYM(J).EQ.'PL') THEN 
                  AC(I,J)=APP
                  AC(J,I)=APP
               ELSE IF (ZSYM(J).EQ.'MI') THEN 
                  AC(I,J)=APM
                  AC(J,I)=APM
               ENDIF
            ELSE IF (ZSYM(I).EQ.'MI') THEN
               IF (ZSYM(J).EQ.'MI') THEN 
                  AC(I,J)=AMM
                  AC(J,I)=AMM
               ELSE IF (ZSYM(J).EQ.'PL') THEN 
                  AC(I,J)=APM
                  AC(J,I)=APM
               ENDIF
            ENDIF 
         ENDDO
      ENDDO
C
C  Calculate the energy.
C
      CALL WENERGY(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV)

      IF (.NOT.GTEST) RETURN
C
C  Gradient.
C
      CALL WGRAD(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,V,XMAT,XMINV)

C     PRINT*,'Numerical derivatives:'
C     DIF=1.0D-4
C     DO J1=1,3*N
C        TEMP1=X(J1)
C        X(J1)=X(J1)+DIF
C        CALL WENERGY(XMU1,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,V1,AC,XMAT,XMINV)
C        X(J1)=X(J1)-2.0D0*DIF
C        CALL WENERGY(XMU2,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,V2,AC,XMAT,XMINV)
C        V(J1)=(V1-V2)/(2.0D0*DIF)
C        WRITE(*,'(I3,F20.10)') J1,V(J1)
C        X(J1)=TEMP1
C     ENDDO

      IF (.NOT.STEST) RETURN
C
C  Hessian.
C
      RETURN
      END
C
C*******************************************************************************
C
C  Energy for the Welch potential
C
      SUBROUTINE WENERGY(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION ESELF, ECC, EREP, ECID, EIDID, POTEL, XMU(3*NATOMS), ALPHA(NATOMS), Q(NATOMS),
     1                 RRD(NATOMS,NATOMS), X(3*NATOMS), XQ(NATOMS), DUMMY, AC(NATOMS,NATOMS),
     2                 RRD3(NATOMS,NATOMS), XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS),
     3                 RRD5(NATOMS,NATOMS), RD, RHOL
      INTEGER I, J, N, K, L

      N=NATOMS
      RHOL=1.0D0/RHO
C
C  Calculate the interparticle distances.
C 
      DO I=1,N
         K=1+3*(I-1)         
         RRD(I,I)=0.0D0
         RRD3(I,I)=0.0D0
         RRD5(I,I)=0.0D0
         DO J=I+1,N
            L=K+3*(J-I)
            RD= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RRD(I,J)=1.0D0/RD
            RRD3(I,J)=RRD(I,J)**3
            RRD5(I,J)=RRD3(I,J)*RRD(I,J)**2
            RRD(J,I)=RRD(I,J)
            RRD3(J,I)=RRD3(I,J)
            RRD5(J,I)=RRD5(I,J)
         ENDDO
      ENDDO
C
C  Calculate the induced dipoles.
C
      CALL DIP(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV)

      ESELF=0.0D0
      ECC=0.0D0
      EREP=0.0D0
      ECID=0.0D0
      EIDID=0.0D0

      DO I=1,N
C
C  Induced dipole self-energy
C
         ESELF=ESELF+(XMU(3*(I-1)+1)**2+XMU(3*(I-1)+2)**2+XMU(3*(I-1)+3)**2)/(2.0D0*ALPHA(I))
         DO J=I+1,N
C
C  Charge-charge.
C
            ECC=ECC+Q(I)*Q(J)*RRD(J,I)
C
C  Exponential repulsion.
C
            DUMMY=(X(3*(J-1)+1)-X(3*(I-1)+1)+XMU(3*(I-1)+1)/XQ(I)-XMU(3*(J-1)+1)/XQ(J))**2
     1           +(X(3*(J-1)+2)-X(3*(I-1)+2)+XMU(3*(I-1)+2)/XQ(I)-XMU(3*(J-1)+2)/XQ(J))**2
     2           +(X(3*(J-1)+3)-X(3*(I-1)+3)+XMU(3*(I-1)+3)/XQ(I)-XMU(3*(J-1)+3)/XQ(J))**2
            DUMMY=SQRT(DUMMY)
            EREP=EREP+AC(J,I)*DEXP(-RHOL*DUMMY)
C
C  Charge-induced dipole term. Include I,J and J,I.
C
            DUMMY=(((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(I-1)+1)+
     1              (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(I-1)+2)+
     2              (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(I-1)+3))*Q(J)
     3           + ((X(3*(I-1)+1)-X(3*(J-1)+1))*XMU(3*(J-1)+1)+
     1              (X(3*(I-1)+2)-X(3*(J-1)+2))*XMU(3*(J-1)+2)+
     2              (X(3*(I-1)+3)-X(3*(J-1)+3))*XMU(3*(J-1)+3))*Q(I))*RRD3(J,I)
            ECID=ECID+DUMMY
C
C  Induced dipole-induced dipole term.
C
            EIDID=EIDID+(XMU(3*(I-1)+1)*XMU(3*(J-1)+1)
     1                  +XMU(3*(I-1)+2)*XMU(3*(J-1)+2)
     2                  +XMU(3*(I-1)+3)*XMU(3*(J-1)+3))*RRD3(J,I)
            DUMMY=((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(I-1)+1)+
     1             (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(I-1)+2)+
     2             (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(I-1)+3))
     3           *((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(J-1)+1)+
     1             (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(J-1)+2)+
     2             (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(J-1)+3))
            EIDID=EIDID-3.0D0*DUMMY*RRD5(J,I)
         ENDDO
      ENDDO
      POTEL=ESELF+ECC+EREP+ECID+EIDID
C     WRITE(*,'(A,F20.10)') 'POTEL=',POTEL
C     WRITE(*,'(A,F20.10)') 'ESELF=',ESELF
C     WRITE(*,'(A,F20.10)') 'ECC=',ECC
C     WRITE(*,'(A,F20.10)') 'EREP=',EREP
C     WRITE(*,'(A,F20.10)') 'ECID=',ECID
C     WRITE(*,'(A,F20.10)') 'EIDID=',EIDID

      RETURN
      END
C
C***************************************************************
C
C  Calculate the induced dipoles by matrix inversion.
C
      SUBROUTINE DIP(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMSAVE,XMINV)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION ALPHA(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS), DUMMY1, XVEC(3*NATOMS),
     2                 DUMMY2, DUMMY3, DUM, D, X(3*NATOMS), Q(NATOMS), XMSAVE(3*NATOMS,3*NATOMS),
     3                 RRD5(NATOMS,NATOMS)
      INTEGER INDX(3*NATOMS), J1, J2, J3, J4
      LOGICAL SFLAG
C
C  Set up the matrix and vector.
C
      DO J1=1,NATOMS
         J3=3*(J1-1)
         XMAT(J3+1,J3+1)=1.0D0
         XMAT(J3+1,J3+2)=0.0D0
         XMAT(J3+1,J3+3)=0.0D0
         XMAT(J3+2,J3+1)=0.0D0
         XMAT(J3+2,J3+2)=1.0D0
         XMAT(J3+2,J3+3)=0.0D0
         XMAT(J3+3,J3+1)=0.0D0
         XMAT(J3+3,J3+2)=0.0D0
         XMAT(J3+3,J3+3)=1.0D0
         DO J2=J1+1,NATOMS
            J4=3*(J2-1)
            DUMMY1=RRD3(J2,J1)
            DUMMY2=RRD5(J2,J1)
            XMAT(J4+1,J3+1)=3.0D0*(X(J4+1)-X(J3+1))**2 * DUMMY2 - DUMMY1
            XMAT(J4+1,J3+2)=3.0D0*(X(J4+1)-X(J3+1))*(X(J4+2)-X(J3+2))*DUMMY2
            XMAT(J4+1,J3+3)=3.0D0*(X(J4+1)-X(J3+1))*(X(J4+3)-X(J3+3))*DUMMY2
            XMAT(J4+2,J3+1)=XMAT(J4+1,J3+2)
            XMAT(J4+2,J3+2)=3.0D0*(X(J4+2)-X(J3+2))**2 * DUMMY2 - DUMMY1
            XMAT(J4+2,J3+3)=3.0D0*(X(J4+2)-X(J3+2))*(X(J4+3)-X(J3+3))*DUMMY2
            XMAT(J4+3,J3+1)=XMAT(J4+1,J3+3)
            XMAT(J4+3,J3+2)=XMAT(J4+2,J3+3)
            XMAT(J4+3,J3+3)=3.0D0*(X(J4+3)-X(J3+3))**2 * DUMMY2 - DUMMY1

            XMAT(J3+1,J4+1)=-XMAT(J4+1,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+1,J4+2)=-XMAT(J4+1,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+1,J4+3)=-XMAT(J4+1,J3+3)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+1)=-XMAT(J4+2,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+2)=-XMAT(J4+2,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+3)=-XMAT(J4+2,J3+3)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+1)=-XMAT(J4+3,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+2)=-XMAT(J4+3,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+3)=-XMAT(J4+3,J3+3)*ALPHA(J1)*0.5D0

            XMAT(J4+1,J3+1)=-XMAT(J4+1,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+1,J3+2)=-XMAT(J4+1,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+1,J3+3)=-XMAT(J4+1,J3+3)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+1)=-XMAT(J4+2,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+2)=-XMAT(J4+2,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+3)=-XMAT(J4+2,J3+3)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+1)=-XMAT(J4+3,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+2)=-XMAT(J4+3,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+3)=-XMAT(J4+3,J3+3)*ALPHA(J2)*0.5D0
         ENDDO
      ENDDO
      
C     PRINT*,'XMAT:'
C     DO J1=1,3*NATOMS
C        WRITE(*,'(6F12.8)') (XMAT(J1,J2),J2=1,3*NATOMS)
C     ENDDO
   
      DO J1=1,NATOMS
         DUMMY1=0.0D0
         DUMMY2=0.0D0
         DUMMY3=0.0D0
         DO J2=1,NATOMS
            DUM=RRD3(J2,J1)*Q(J2)
            DUMMY1=DUMMY1-(X(3*(J1-1)+1)-X(3*(J2-1)+1))*DUM
            DUMMY2=DUMMY2-(X(3*(J1-1)+2)-X(3*(J2-1)+2))*DUM
            DUMMY3=DUMMY3-(X(3*(J1-1)+3)-X(3*(J2-1)+3))*DUM
         ENDDO
         XVEC(3*(J1-1)+1)=-DUMMY1*ALPHA(J1)*0.5D0
         XVEC(3*(J1-1)+2)=-DUMMY2*ALPHA(J1)*0.5D0
         XVEC(3*(J1-1)+3)=-DUMMY3*ALPHA(J1)*0.5D0
      ENDDO

C     PRINT*,'XVEC:'
C     DO J1=1,3*NATOMS
C        WRITE(*,'(I3,F12.8)') J1,XVEC(J1)
C     ENDDO

C
C  Initialize XMINV to the unit matrix and save XMAT.
C
      DO J1=1,3*NATOMS
         XMINV(J1,J1)=1.0D0
         XMSAVE(J1,J1)=XMAT(J1,J1)
         DO J2=J1+1,3*NATOMS
            XMINV(J2,J1)=0.0D0
            XMINV(J1,J2)=0.0D0
            XMSAVE(J2,J1)=XMAT(J2,J1)
            XMSAVE(J1,J2)=XMAT(J1,J2)
         ENDDO
      ENDDO
      
      CALL LUDCMP(XMAT,3*NATOMS,3*NATOMS,INDX,D,SFLAG)
      IF (SFLAG) THEN
         PRINT*,'Singular matrix'
      ENDIF

      DO J1=1,3*NATOMS
         CALL LUBKSB(XMAT,3*NATOMS,3*NATOMS,INDX,XMINV(1,J1))
      ENDDO

C     PRINT*,'XMINV:'
C     DO J1=1,3*NATOMS
C        WRITE(*,'(6E20.10)') (XMINV(J1,J2),J2=1,3*NATOMS)
C     ENDDO

      DO J1=1,3*NATOMS
         DUMMY1=0.0D0
         DO J2=1,3*NATOMS
            DUMMY1=DUMMY1+XMINV(J1,J2)*XVEC(J2)
         ENDDO
         XMU(J1)=DUMMY1
         PRINT*,'J1,XMU=',J1,XMU(J1)
      ENDDO

      RETURN
      END
C
C***************************************************************
C
C  Calculate the analytic derivatives of the induced dipoles.
C
      SUBROUTINE DIPGRAD(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV,XMUGRAD)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION ALPHA(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS), DUMMY, 
     2                 DUM, X(3*NATOMS), Q(NATOMS), VEC3(3*NATOMS), RRD5(NATOMS,NATOMS),
     3                 XMUGRAD(3*NATOMS,3*NATOMS), VEC1(3*NATOMS), VEC2(3*NATOMS)
      INTEGER J1, J2, J3, J4, J5, J6, J7
      
C
C  Set up the derivative matrix and vector for atom J5 component J6.
C
C     PRINT*,'Analytic derivatives of mu'
      DO J5=1,NATOMS           ! J5 = i (can be A or B)
         J4=3*(J5-1)
         DO J6=1,3             ! J6 = gamma
            DO J1=1,3*NATOMS
               VEC2(J1)=0.0D0
            ENDDO
            DO J7=1,3          ! first index of T tensor = alpha
               DO J1=1,NATOMS  ! atom B
                  J3=3*(J1-1)
                  DUM=-15.0D0*RRD5(J1,J5)*ALPHA(J5)*(X(J4+J7)-X(J3+J7))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2
                  DO J2=1,3
                     VEC1(J3+J2)=(X(J4+J2)-X(J3+J2))*DUM
                  ENDDO
                  DUM=RRD5(J1,J5)*ALPHA(J5)
                  VEC1(J3+J7)=VEC1(J3+J7)+3.0D0*(X(J4+J6)-X(J3+J6))*DUM
                  VEC1(J3+J6)=VEC1(J3+J6)+3.0D0*(X(J4+J7)-X(J3+J7))*DUM
               ENDDO
               DUMMY=0.0D0
               DO J1=1,3*NATOMS
                  DUMMY=DUMMY+VEC1(J1)*XMU(J1)
               ENDDO
               VEC2(J4+J7)=DUMMY
            ENDDO
            DO J1=1,NATOMS  ! atom B
               J3=3*(J1-1)
               DUM=RRD5(J1,J5)
               VEC1(J3+1)=(X(J4+1)-X(J3+1))*DUM
               VEC1(J3+2)=(X(J4+2)-X(J3+2))*DUM
               VEC1(J3+3)=(X(J4+3)-X(J3+3))*DUM
            ENDDO
            DUMMY=0.0D0
            DO J1=1,3*NATOMS
               DUMMY=DUMMY+VEC1(J1)*XMU(J1)
            ENDDO
            VEC2(J4+J6)=VEC2(J4+J6)+DUMMY*3.0D0*ALPHA(J5)
   
            DO J7=1,3          ! first index of T tensor = alpha
               DO J1=1,NATOMS  ! atom A this time
                  J3=3*(J1-1)
                  DUM=RRD5(J1,J5)*ALPHA(J1)
                  DUMMY=(X(J4+1)-X(J3+1))*XMU(J4+1)+(X(J4+2)-X(J3+2))*XMU(J4+2)+(X(J4+3)-X(J3+3))*XMU(J4+3)
                  VEC2(J3+J7)=VEC2(J3+J7)+3.0D0*DUM*((X(J4+J6)-X(J3+J6))*XMU(J4+J7)
     1                                              +(X(J4+J7)-X(J3+J7))*(XMU(J4+J6)
     2                                  -5.0D0*DUMMY*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2))
               ENDDO
            ENDDO
            DO J1=1,NATOMS  ! atom A this time
               J3=3*(J1-1)
               DUM=RRD5(J1,J5)*ALPHA(J1)
               DO J2=1,3    ! beta index of T(Ai)(2) where B or i=J1
                  DUMMY=3.0D0*(X(J4+J2)-X(J3+J2))*DUM
                  VEC2(J3+J6)=VEC2(J3+J6)+DUMMY*XMU(J4+J2)
               ENDDO
            ENDDO
C
C  Now VEC2 should contain (d M/d X) mu apart from the factor of -1/2.
C
            DO J1=1,NATOMS   ! second atom B
               J3=3*(J1-1)
               DUM=RRD3(J1,J5)*ALPHA(J1)*Q(J5)
               DO J2=1,3     ! alpha index of T(iB)(1)
                  VEC3(J3+J2)=-3.0D0*(X(J4+J2)-X(J3+J2))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2*DUM
               ENDDO
               VEC3(J3+J6)=VEC3(J3+J6)+DUM
            ENDDO

            DO J1=1,NATOMS   ! first atom A
               J3=3*(J1-1)
               DUM=RRD3(J1,J5)*ALPHA(J5)*Q(J1)
               DO J2=1,3     ! alpha index of T(Ai)(1)
                  VEC3(J4+J2)= VEC3(J4+J2)+3.0D0*(X(J4+J2)-X(J3+J2))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2*DUM
               ENDDO
               VEC3(J4+J6)=VEC3(J4+J6)-DUM
            ENDDO
C
C  Now VEC3 - VEC2 should contain (d Y/d X) - (d M/d X) mu apart from a factor of -1/2.
C
            DO J1=1,3*NATOMS
               DUMMY=0.0D0
               DO J2=1,3*NATOMS
                  DUMMY=DUMMY-0.5D0*XMINV(J1,J2)*(VEC3(J2)-VEC2(J2))
               ENDDO
               XMUGRAD(J1,J4+J6)=DUMMY
C              WRITE(*,'(2I3,F20.10)') J1,J4+J6,XMUGRAD(J1,J4+J6)
            ENDDO

         ENDDO
      ENDDO

      RETURN
      END

C
C*******************************************************************************
C
C  Analytic gradient for the Welch potential
C
      SUBROUTINE WGRAD(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,GRAD,XMAT,XMINV)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION VSELF(3*NATOMS), VCC(3*NATOMS), VREP(3*NATOMS), VCID(3*NATOMS), VIDID(3*NATOMS), 
     1                 GRAD(3*NATOMS), XMU(3*NATOMS), ALPHA(NATOMS), Q(NATOMS),
     2                 RRD(NATOMS,NATOMS), X(3*NATOMS), XQ(NATOMS), DUMMY, AC(NATOMS,NATOMS),
     3                 XMUGRAD(3*NATOMS,3*NATOMS), XMAT(3*NATOMS,3*NATOMS),
     4                 XMINV(3*NATOMS,3*NATOMS), DUM, DUM2, VD(3), RDUM, RRD3(NATOMS,NATOMS), RRD5(NATOMS,NATOMS),
     5                 T1, T2, T3, RHOL
      INTEGER N, J1, J2, J3, J4, J5, J6, J7

      N=NATOMS
      RHOL=1.0D0/RHO
C
C  Calculate induced dipole derivatives.
C
      CALL DIPGRAD(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV,XMUGRAD)
C
C  Charge-charge.
C
      DO J1=1,N
         DO J2=1,3
            DUMMY=0.0D0
            J4=3*(J1-1)+J2
            DO J3=1,N
               DUMMY=DUMMY-Q(J3)*RRD3(J3,J1)*(X(J4)-X(3*(J3-1)+J2))
            ENDDO
            VCC(J4)=Q(J1)*DUMMY
         ENDDO
      ENDDO
C
C  Self-energy.
C
      DO J1=1,N       ! this is i
         DO J3=1,3    ! this is alpha
            J2=3*(J1-1)+J3
            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J5+1)*XMUGRAD(J5+1,J2)
     1                     +XMU(J5+2)*XMUGRAD(J5+2,J2)
     2                     +XMU(J5+3)*XMUGRAD(J5+3,J2))/ALPHA(J4)
            ENDDO
            VSELF(J2)=DUMMY
         ENDDO
      ENDDO
C
C  Charge-induced dipole
C
      DO J1=1,N         ! atom i
         J2=3*(J1-1)
         DO J3=1,3      ! index gamma
            DUMMY=0.0D0
            DO J4=1,N   ! A
               DUMMY=DUMMY+XMU(3*(J4-1)+J3)*RRD3(J4,J1)
            ENDDO
            VCID(J2+J3)=Q(J1)*DUMMY
            DUMMY=0.0D0
            DO J4=1,N   ! B
               DUMMY=DUMMY-Q(J4)*RRD3(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)+DUMMY*XMU(J2+J3)
            DUMMY=0.0D0
            DO J4=1,N   ! B
               J5=3*(J4-1)
               DUMMY=DUMMY+Q(J4)*(X(J2+J3)-X(3*(J4-1)+J3))*
     1               ((X(J5+1)-X(J2+1))*XMU(J2+1)
     2               +(X(J5+2)-X(J2+2))*XMU(J2+2)
     3               +(X(J5+3)-X(J2+3))*XMU(J2+3))*RRD5(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)-3.0D0*DUMMY
            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(X(J2+J3)-X(J5+J3))*
     1               ((X(J2+1)-X(J5+1))*XMU(J5+1)
     2               +(X(J2+2)-X(J5+2))*XMU(J5+2)
     3               +(X(J2+3)-X(J5+3))*XMU(J5+3))*RRD5(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)-3.0D0*DUMMY*Q(J1)
            DUMMY=0.0D0
            DO J4=1,N    ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N ! B
                  J7=3*(J5-1)
                  DUMMY=DUMMY+Q(J5)*((X(J7+1)-X(J6+1))*T1
     1                              +(X(J7+2)-X(J6+2))*T2
     2                              +(X(J7+3)-X(J6+3))*T3)*RRD3(J5,J4)
               ENDDO
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)+DUMMY
         ENDDO
      ENDDO
C
C  Induced dipole-induced dipole.
C 
      DO J1=1,N
         J2=3*(J1-1)
         DO J3=1,3
            DUMMY=0.0D0
            DO J4=1,N
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N
                  J7=3*(J5-1)
                  DUMMY=DUMMY+(XMU(J7+1)*T1+XMU(J7+2)*T2+XMU(J7+3)*T3)*RRD3(J5,J4)
               ENDDO
            ENDDO
            VIDID(J2+J3)=DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J2+1)*XMU(J5+1)+XMU(J2+2)*XMU(J5+2)+XMU(J2+3)*XMU(J5+3))
     1                    *(X(J2+J3)-X(J5+J3))*RRD5(J1,J4)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N     ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N  ! B
                  J7=3*(J5-1)
                  DUMMY=DUMMY+(T1*(X(J7+1)-X(J6+1))
     1                        +T2*(X(J7+2)-X(J6+2))
     2                        +T3*(X(J7+3)-X(J6+3)))
     3                       *(XMU(J7+1)*(X(J7+1)-X(J6+1))
     4                        +XMU(J7+2)*(X(J7+2)-X(J6+2))
     5                        +XMU(J7+3)*(X(J7+3)-X(J6+3)))*RRD5(J5,J4)
               ENDDO
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+XMU(J5+J3)*(XMU(J2+1)*(X(J2+1)-X(J5+1))
     1                                +XMU(J2+2)*(X(J2+2)-X(J5+2))
     2                                +XMU(J2+3)*(X(J2+3)-X(J5+3)))*RRD5(J4,J1)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J5+1)*(X(J2+1)-X(J5+1))
     1                     +XMU(J5+2)*(X(J2+2)-X(J5+2))
     2                     +XMU(J5+3)*(X(J2+3)-X(J5+3)))*RRD5(J4,J1)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY*XMU(J2+J3)

            DUMMY=0.0D0
            DO J4=1,N  ! B
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J2+1)*(X(J2+1)-X(J5+1))
     1                     +XMU(J2+2)*(X(J2+2)-X(J5+2))
     2                     +XMU(J2+3)*(X(J2+3)-X(J5+3)))
     3                    *(XMU(J5+1)*(X(J2+1)-X(J5+1))
     4                     +XMU(J5+2)*(X(J2+2)-X(J5+2))
     5                     +XMU(J5+3)*(X(J2+3)-X(J5+3)))
     6                    *(X(J2+J3)-X(J5+J3))*RRD5(J4,J1)*RRD(J4,J1)**2
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)+15.0D0*DUMMY
         ENDDO
      ENDDO
C
C  Exponential repulsion.
C
      DO J1=1,N
         J2=3*(J1-1)
         DO J3=1,3
            DUMMY=0.0D0
            DO J4=1,N     ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=J4+1,N  ! B
                  J7=3*(J5-1)
                  DUM=(X(J7+1)-X(J6+1)+XMU(J6+1)/XQ(J4)-XMU(J7+1)/XQ(J5))**2
     1               +(X(J7+2)-X(J6+2)+XMU(J6+2)/XQ(J4)-XMU(J7+2)/XQ(J5))**2
     2               +(X(J7+3)-X(J6+3)+XMU(J6+3)/XQ(J4)-XMU(J7+3)/XQ(J5))**2
                  DUM=SQRT(DUM)
                  RDUM=1.0D0/DUM
                  VD(1)=T1/XQ(J4)-XMUGRAD(J7+1,J2+J3)/XQ(J5)
                  VD(2)=T2/XQ(J4)-XMUGRAD(J7+2,J2+J3)/XQ(J5)
                  VD(3)=T3/XQ(J4)-XMUGRAD(J7+3,J2+J3)/XQ(J5)
                  IF (J1.EQ.J4) VD(J3)=VD(J3)-1.0D0
                  IF (J1.EQ.J5) VD(J3)=VD(J3)+1.0D0
                  DUM2=VD(1)*(X(J7+1)-X(J6+1)+XMU(J6+1)/XQ(J4)-XMU(J7+1)/XQ(J5))
     1                +VD(2)*(X(J7+2)-X(J6+2)+XMU(J6+2)/XQ(J4)-XMU(J7+2)/XQ(J5))
     2                +VD(3)*(X(J7+3)-X(J6+3)+XMU(J6+3)/XQ(J4)-XMU(J7+3)/XQ(J5))
                  DUMMY=DUMMY+AC(J5,J4)*DEXP(-RHOL*DUM)*DUM2*RDUM
               ENDDO
            ENDDO
            VREP(J2+J3)=-DUMMY*RHOL

         ENDDO
      ENDDO

      DO J1=1,3*N
         GRAD(J1)=VCC(J1)+VSELF(J1)+VREP(J1)+VCID(J1)+VIDID(J1)
      ENDDO
C     PRINT*,'Analytic gradient:'
C     WRITE(*,'(I3,F20.10)') (J1,GRAD(J1),J1=1,3*NATOMS)

      RETURN
      END

C
C***************************************************************
C
C  Matrix inversion subroutines.
C

      SUBROUTINE LUDCMP(A,N,NP,INDX,D,SFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=1000,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(NP),VV(NMAX)
      LOGICAL SFLAG
      IF (N.GT.NMAX) THEN
         PRINT*,'*** WARNING - increase NMAX in ludcmp!!'
      ENDIF
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) THEN
           PRINT*,'Singular matrix.'
           SFLAG=.TRUE.
           RETURN
        ENDIF
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.)A(J,J)=TINY
        IF(J.NE.N)THEN
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      RETURN
      END
 
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(NP),B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            ENDDO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         ENDIF
         B(I)=SUM
      ENDDO

      DO I=N,1,-1
         SUM=B(I)
         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            ENDDO
         ENDIF
         B(I)=SUM/A(I,I)
      ENDDO

      RETURN
      END
