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
      SUBROUTINE AMB(XVEC,GRAD,EREAL,GRADT)
      USE commons
      USE modamber  
      IMPLICIT NONE
      DOUBLE PRECISION XVEC(3*NATOMS), GRAD(3*NATOMS), EREAL
      DOUBLE PRECISION vixen1, vixen2, vixen3, vixen4, vixen5, vixen6, e1, e2, e3 
      LOGICAL GRADT
      INTEGER J1, J2
      SAVE

      IF (ang.EQ.0) THEN
        count=1
        IF (MGBWATER) CALL GETALPHA
        CALL AMBERMASS
        CALL AMBERS
      END IF
      ambercount= ambercount+1

      IF (count.EQ.1) THEN
        DO J1=1,NATOMS
          XVEC(3*(J1-1)+1)=x(J1)
          XVEC(3*(J1-1)+2)=y(J1)
          XVEC(3*(J1-1)+3)=z(J1)
        ENDDO
      ELSE
        DO J1=1,NATOMS
           x(J1)=XVEC(3*(J1-1)+1)
           y(J1)=XVEC(3*(J1-1)+2)
           z(J1)=XVEC(3*(J1-1)+3)
        ENDDO
      END IF
      count=0

      CALL amberenergy

      EREAL=totenergy
C     WRITE (*,'(F20.10)') EREAL

      IF (GRADT) THEN
         ambergcount= ambergcount+1
         CALL AMBG
         IF (FIX) THEN
C            WRITE (*,'(A)') 'Fixing backbone atoms and side-chain chiral centres'
            DO J1=1,NATOMS
               IF (label(J1).EQ.'d') THEN
                  dEbydx(J1)=0.0D0
                  dEbydy(J1)=0.0D0
                  dEbydz(J1)=0.0D0
                  DO J2=1,NATOMS
                     IF (bonds(J1,J2).EQ.1) THEN
                        dEbydx(J2)=0.0D0
                        dEbydy(J2)=0.0D0
                        dEbydz(J2)=0.0D0
                     END IF
                  END DO
               END IF
            END DO
         END IF
         arms=0.0D0
         DO J1=1,NATOMS
            GRAD(3*(J1-1)+1)=dEbydx(J1)
            GRAD(3*(J1-1)+2)=dEbydy(J1)
            GRAD(3*(J1-1)+3)=dEbydz(J1)
         ENDDO
      ENDIF


C      CALL numergrad

C      CALL secondderivs

C      CALL numersecderiv


      RETURN
      END

      SUBROUTINE one4
      USE modamber
      IMPLICIT NONE

      DO a=1,t
        one_four(da1(a),da4(a))=1
      ENDDO

C     PRINT *,"One-four relationships"
C     DO a=1,atoms
C       DO b=a+1,atoms
C         IF (one_four(a,b).EQ.1)  PRINT *,a,"-",b
C       ENDDO
C     ENDDO

      RETURN
      END
C
C  Setup stuff only needs doing once.
C
      SUBROUTINE AMBERS
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER            marvin
      INTEGER            chiratom(4)
      INTEGER            canine
      COMMON /CHIR/      canine

      ambercount=0
      ambergcount=0
      NDIHEDRALS=0
C 
C Create entry in angarray
C 
      DO b=1,atoms
         DO c=1,atoms
            DO d=b,atoms
               IF (b.NE.c .AND. c.NE.d .AND. b.NE.d) THEN
                  IF (.NOT.(bonds(b,c).NE.1 .OR. bonds(c,d).NE.1)) THEN
                     ang=ang+1
                     aa1(ang)=b
                     aa2(ang)=c
                     aa3(ang)=d
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO i=1,atoms
         DO j=1,atoms
            DO k=1,atoms
               DO l=i,atoms
                  IF (bonds(i,j).NE.1 .OR. bonds(j,i).NE.1) GOTO 113
                  IF (bonds(j,k).NE.1 .OR. bonds(k,j).NE.1) GOTO 113
                  IF (bonds(k,l).NE.1 .OR. bonds(l,k).NE.1) GOTO 113
                  IF (i.EQ.j .OR. i.EQ.k .OR. i.EQ.l .OR. j.EQ.k .OR. j.EQ.l .OR. k.EQ.l) GOTO 113
  
                  colin=0
                  ax=x(j)-x(i)
                  ay=y(j)-y(i)
                  az=z(j)-z(i)
                  bx=x(k)-x(j)
                  by=y(k)-y(j)
                  bz=z(k)-z(j)
                  cx=x(j)-x(l)
                  cy=y(j)-y(l)
                  cz=z(j)-z(l)
                  dx=x(l)-x(k)
                  dy=y(l)-y(k)
                  dz=z(l)-z(k)
C 
C Checking for colinearity,connectivity and coincidence
C 
                  IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) 
     1               .AND. ay*bz.GT.(by*az-TINY) 
     2               .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
                  IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) 
     1               .AND. by*dz.GT.(dy*bz-TINY) 
     2               .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2
                  IF (colin.EQ.1) THEN
                    PRINT *,'Three sites',i,j,k,'are colinear'
                    STOP
                  ELSE IF (colin.EQ.2) THEN
                    PRINT *,'Three sites',j,k,l,'are colinear'
                    STOP
                  ELSE
C 
C Create entry in torsarray
C
                    t=t+1
                    da1(t)=i
                    da2(t)=j
                    da3(t)=k
                    da4(t)=l
                  ENDIF
113               CONTINUE
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL one4

      PRINT *,' Torsions assigned'

      DO a=1,t
         b=type(da1(a))
         c=type(da2(a))
         d=type(da3(a))        
         e=type(da4(a))
         match=0
         DO f=1,60
            IF (gentorsparams(f,1).LT.0.9) GOTO 5
            g=gentorsparams(f,1)
            h=gentorsparams(f,2)
            IF (c.EQ.g .AND. d.EQ.h) THEN
               match=1
               did(a)=gentorsparams(f,3)
               dvn(a)=gentorsparams(f,4)
               ddelta(a)=gentorsparams(f,5)
               dn(a)=gentorsparams(f,6)
               dvn2(a)=0.0
               dvn3(a)=0.0
               ddelta2(a)=0.0
               ddelta3(a)=0.0
               dn2(a)=0.0
               dn3(a)=0.0
            ELSE IF (d.EQ.g .AND. c.EQ.h) THEN
               match=1
               did(a)=gentorsparams(f,3)
               dvn(a)=gentorsparams(f,4)
               ddelta(a)=gentorsparams(f,5)
               dn(a)=gentorsparams(f,6)
               dvn2(a)=0.0
               dvn3(a)=0.0
               ddelta2(a)=0.0
               ddelta3(a)=0.0
               dn2(a)=0.0
               dn3(a)=0.0
            ENDIF
         ENDDO
5        CONTINUE

         DO f=1,50
            IF (spectorsparams(f,1).LT.0.9) GOTO 6
            g=spectorsparams(f,1)
            h=spectorsparams(f,2)
            i=spectorsparams(f,3)
            j=spectorsparams(f,4)
            IF (b.EQ.g .AND. c.EQ.h .AND. d.EQ.i .AND. e.EQ.j) THEN
               match=1
               did(a)=spectorsparams(f,5)
               dvn(a)=spectorsparams(f,6)
               ddelta(a)=spectorsparams(f,7)
               dn(a)=spectorsparams(f,8)
               dvn2(a)=spectorsparams(f,9)
               ddelta2(a)=spectorsparams(f,10)
               dn2(a)=spectorsparams(f,11)
               dvn3(a)=spectorsparams(f,12)
               ddelta3(a)=spectorsparams(f,13)
               dn3(a)=spectorsparams(f,14)
            ELSE IF (e.EQ.g .AND. d.EQ.h .AND. c.EQ.i .AND. b.EQ.j) THEN
               match=1
               did(a)=spectorsparams(f,5)
               dvn(a)=spectorsparams(f,6)
               ddelta(a)=spectorsparams(f,7)
               dn(a)=spectorsparams(f,8)
               dvn2(a)=spectorsparams(f,9)
               ddelta2(a)=spectorsparams(f,10)
               dn2(a)=spectorsparams(f,11)
               dvn3(a)=spectorsparams(f,12)
               ddelta3(a)=spectorsparams(f,13)
               dn3(a)=spectorsparams(f,14)
            ENDIF
         ENDDO
6        CONTINUE
      ENDDO

      IF (match.EQ.0) THEN
         did(a)=1
         dvn(a)=0.0
         dvn2(a)=0.0
         dvn3(a)=0.0
      ENDIF

      DO i=1,atoms
       DO j=i+1,atoms
        DO k=j+1,atoms
         DO l=k+1,atoms
          IF (.NOT.(i.EQ.j .OR. i.EQ.k .OR. i.EQ.l .OR. j.EQ.k .OR. j.EQ.l .OR. k.EQ.l)) THEN
C 
C This IF construct puts the atoms in the improper in the correct order 1-2-3-4
C 
          IF (bonds(i,k).EQ.1 .AND. bonds(j,k).EQ.1 .AND. bonds(k,l).EQ.1) THEN
           imp=imp+1
           IF ((type(i) .LE. type(j)) .AND. (type(j) .LE. type(l))) THEN 
            ia1(imp)=i
            ia2(imp)=j
            ia3(imp)=k
            ia4(imp)=l
           ELSE IF ((type(i) .LE. type(l)) .AND. (type(l) .LT. type(j))) THEN
            ia1(imp)=i
            ia2(imp)=l
            ia3(imp)=k
            ia4(imp)=j
           ELSE IF ((type(j) .LT. type(i)) .AND. (type(i) .LE. type(l))) THEN
            ia1(imp)=j
            ia2(imp)=i
            ia3(imp)=k
            ia4(imp)=l
           ELSE IF ((type(j) .LE. type(l)) .AND. (type(l) .LT. type(i))) THEN
            ia1(imp)=j
            ia2(imp)=l
            ia3(imp)=k
            ia4(imp)=i
           ELSE IF ((type(l) .LT. type(i)) .AND. (type(i) .LE. type(j))) THEN
            ia1(imp)=l
            ia2(imp)=i
            ia3(imp)=k
            ia4(imp)=j
           ELSE IF ((type(l) .LT. type(j)) .AND. (type(j) .LT. type(i))) THEN
            ia1(imp)=l
            ia2(imp)=j
            ia3(imp)=k
            ia4(imp)=i
           ENDIF
            
          ELSE IF (bonds(i,j).EQ.1 .AND. bonds(i,k).EQ.1 .AND. bonds(i,l).EQ.1) THEN
           imp=imp+1
           IF ((type(j) .LE. type(k)) .AND. (type(k) .LE. type(l))) THEN 
            ia1(imp)=j
            ia2(imp)=k
            ia3(imp)=i
            ia4(imp)=l
           ELSE IF ((type(j) .LE. type(l)) .AND. (type(l) .LT. type(k))) THEN
            ia1(imp)=j
            ia2(imp)=l
            ia3(imp)=i
            ia4(imp)=k
           ELSE IF ((type(k) .LT. type(j)) .AND. (type(j) .LE. type(l))) THEN
            ia1(imp)=k
            ia2(imp)=j
            ia3(imp)=i
            ia4(imp)=l
           ELSE IF ((type(k) .LE. type(l)) .AND. (type(l) .LT. type(j))) THEN
            ia1(imp)=k
            ia2(imp)=l
            ia3(imp)=i
            ia4(imp)=j
           ELSE IF ((type(l) .LT. type(j)) .AND. (type(j) .LE. type(k))) THEN
            ia1(imp)=l
            ia2(imp)=j
            ia3(imp)=i
            ia4(imp)=k
           ELSE IF ((type(l) .LT. type(k)) .AND. (type(k) .LT. type(j))) THEN
            ia1(imp)=l
            ia2(imp)=k
            ia3(imp)=i
            ia4(imp)=j
           ENDIF

          ELSE IF (bonds(i,j).EQ.1 .AND. bonds(j,k).EQ.1 .AND. bonds(j,l).EQ.1) THEN
           imp=imp+1
           IF ((type(i) .LE. type(k)) .AND. (type(k) .LE. type(l))) THEN 
            ia1(imp)=i
            ia2(imp)=k
            ia3(imp)=j
            ia4(imp)=l
           ELSE IF ((type(i) .LE. type(l)) .AND. (type(l) .LT. type(k))) THEN
            ia1(imp)=i
            ia2(imp)=l
            ia3(imp)=j
            ia4(imp)=k
           ELSE IF ((type(k) .LT. type(i)) .AND. (type(i) .LE. type(l))) THEN
            ia1(imp)=k
            ia2(imp)=i
            ia3(imp)=j
            ia4(imp)=l
           ELSE IF ((type(k) .LE. type(l)) .AND. (type(l) .LT. type(i))) THEN
            ia1(imp)=k
            ia2(imp)=l
            ia3(imp)=j
            ia4(imp)=i
           ELSE IF ((type(l) .LT. type(i)) .AND. (type(i) .LE. type(k))) THEN
            ia1(imp)=l
            ia2(imp)=i
            ia3(imp)=j
            ia4(imp)=k
           ELSE IF ((type(l) .LT. type(k)) .AND. (type(k) .LT. type(i))) THEN
            ia1(imp)=l
            ia2(imp)=k
            ia3(imp)=j
            ia4(imp)=i
           ENDIF

          ELSE IF (bonds(i,l).EQ.1 .AND. bonds(j,l).EQ.1 .AND. bonds(k,l).EQ.1) THEN
           imp=imp+1
           IF ((type(i) .LE. type(j)) .AND. (type(j) .LE. type(k))) THEN 
            ia1(imp)=i
            ia2(imp)=j
            ia3(imp)=l
            ia4(imp)=k
           ELSE IF ((type(i) .LE. type(k)) .AND. (type(k) .LT. type(j))) THEN
            ia1(imp)=i
            ia2(imp)=k
            ia3(imp)=l
            ia4(imp)=j
           ELSE IF ((type(j) .LT. type(i)) .AND. (type(i) .LE. type(k))) THEN
            ia1(imp)=j
            ia2(imp)=i
            ia3(imp)=l
            ia4(imp)=k
           ELSE IF ((type(j) .LE. type(k)) .AND. (type(k) .LT. type(i))) THEN
            ia1(imp)=j
            ia2(imp)=k
            ia3(imp)=l
            ia4(imp)=i
           ELSE IF ((type(k) .LT. type(i)) .AND. (type(i) .LE. type(j))) THEN
            ia1(imp)=k
            ia2(imp)=i
            ia3(imp)=l
            ia4(imp)=j
           ELSE IF ((type(k) .LT. type(j)) .AND. (type(j) .LT. type(i))) THEN
            ia1(imp)=k
            ia2(imp)=j
            ia3(imp)=l
            ia4(imp)=i
           ENDIF
           ENDIF
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      PRINT *,' Impropers assigned'

C
C Assign parameters to each improper
C
      DO a=1,imp
         b=type(ia1(a))
         c=type(ia2(a))        
         d=type(ia3(a))
         e=type(ia4(a))
         DO f=1,15
           g=genimpparams(f,1)
           h=genimpparams(f,2)
           IF (d.EQ.g .AND. (e.EQ.h .OR. c.EQ.h .OR. b.EQ.h)) THEN
             ivn(a)=genimpparams(f,3)
             idelta(a)=genimpparams(f,4)
             in1(a)=genimpparams(f,5)
           ENDIF
         ENDDO
         DO f=1,4
           g=midimpparams(f,1)
           h=midimpparams(f,2)
           i=midimpparams(f,3)
           IF (c.EQ.g .AND. d.EQ.h .AND. e.EQ.i) THEN
             ivn(a)=midimpparams(f,4)
             idelta(a)=midimpparams(f,5)
             in1(a)=midimpparams(f,6)
           ELSE IF (b.EQ.g .AND. c.EQ.i .AND. d.EQ.h) THEN
             ivn(a)=midimpparams(f,4)
             idelta(a)=midimpparams(f,5)
             in1(a)=midimpparams(f,6)
           ENDIF
         ENDDO
         DO f=1,15
           g=specimpparams(f,1)
           h=specimpparams(f,2)
           i=specimpparams(f,3)
           j=specimpparams(f,4)
   
           IF (b.EQ.g .AND. c.EQ.h .AND. d.EQ.i .AND. e.EQ.j) THEN
             ivn(a)=specimpparams(f,5)
             idelta(a)=specimpparams(f,6)
             in1(a)=specimpparams(f,7)
           ENDIF
         ENDDO
      ENDDO

      DO a=1,atoms
         DO b=a+1,atoms
            rstar=vdwr(type(a))+vdwr(type(b))
            epsilon=SQRT(vdwe(type(a))*vdwe(type(b)))
            vdwa(a,b)=epsilon*rstar**12
            vdwb(a,b)=2*epsilon*rstar**6  
         ENDDO
      ENDDO

      canine=1
      DO a=1,atoms
        IF (chiral(a).EQ.1) GOTO 111

        IF (typech(a).NE.'CT') GOTO 112
        marvin=1
        DO b=1,atoms
          IF (bonds(a,b).EQ.1 .OR. bonds(b,a).EQ.1) THEN
            chiratom(marvin)=b
            marvin=marvin+1
          END IF
        END DO

        i=chiratom(1)
        j=chiratom(2)
        k=chiratom(3)
        l=chiratom(4)

        IF (typech(i).NE.typech(j) .AND. typech(i).NE.typech(k) .AND. typech(i).NE.typech(l) .AND. 
     1      typech(j).NE.typech(k) .AND. typech(j).NE.typech(l) .AND. typech(k).NE.typech(l)) THEN
          chiral(a)=1
        ELSE
          GOTO 112
        END IF
111     CONTINUE

        chiralarray(canine,1)=a
        chiralarray(canine,2)=i
        chiralarray(canine,3)=j
        chiralarray(canine,4)=k
        chiralarray(canine,5)=l
        canine=canine+1
   
112   CONTINUE
      END DO

      PRINT *,'Set up routine completed'
      PRINT *,'No. of angles - ',ang
      PRINT *,'No. of dihedrals - ',t
      PRINT *,'No. of impropers - ',imp

      canine=canine-1
      DO a=1,canine
        i=chiralarray(a,2)
        j=chiralarray(a,3)
        k=chiralarray(a,4)
        l=chiralarray(a,5)
C        PRINT *,'Atom',chiralarray(a,1),' is chiral and is bonded to atoms',i,' ',j,' ',k,' ',l
      END DO

      DO a=1,atoms
        IF (NDIHEDRALS.EQ.0 .AND. label(a).EQ.'t' .AND. typech(a).EQ.'N3') THEN
          CAP=.FALSE.
          WRITE (*,*) 'No capping groups'
        END IF
        IF (label(a).EQ.'d' .AND. (typech(a).NE.'N3' .AND. typech(a).NE.'C')) THEN
          DO b=a+1,atoms
            IF (label(b).EQ.'d' .AND. bonds(a,b).EQ.1) THEN
              NDIHEDRALS=NDIHEDRALS+1
              DATOM1(NDIHEDRALS)=a
              DATOM2(NDIHEDRALS)=b
C              PRINT *,'Dihedral ',NDIHEDRALS,'  ',a,'-',b
C              PRINT *,'         ',typech(a),' - ',typech(b)
            END IF
          END DO
        END IF
      END DO

      CALL makelist

      RETURN
      END
C
C  Distances
C  Van der Waals and charge terms
C
      SUBROUTINE AMBERD
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER NDUMMY1, NDUMMY3
      DOUBLE PRECISION DUMMY, DUMMY2, vixen1, vixen2, vixen3, vixen4, vixen5, fgb, alpha_ij, vixen6
      DOUBLE PRECISION DP1

      IF (AMCUT) THEN
         DO c=1,rlistcount
            a=rlist(c,1)
            b=rlist(c,2)
            vixen4=q(a)
            DUMMY2=(x(a)-x(b))**2+(y(a)-y(b))**2+(z(a)-z(b))**2
            vixen5=SQRT(DUMMY2)

            NDUMMY1=1-bonds(a,b)
            NDUMMY3=NDUMMY1*(1-one_three(a,b))
            IF (NDUMMY3.NE.0) THEN
               vixen3=FLOAT(one_four(a,b))
               DUMMY=1.0D0/DUMMY2**3
               vixen1=(vdwa(a,b)*DUMMY-vdwb(a,b))*DUMMY - (vdwa(a,b)*DP1-vdwb(a,b))*DP1
               vdwenergy=vdwenergy+vixen1*(1.0D0-vixen3/2.0D0)
               IF (FAKEWATER) THEN
                  vixen2=vixen4*q(b)/(dparam*dielec*vixen5**2) - vixen4*q(b)/(dparam*dielec*QCUTOFF**2)
               ELSE
                  vixen2=vixen4*q(b)/(vixen5*dielec*dparam) - vixen4*q(b)/(dparam*dielec*QCUTOFF)
               END IF
               qenergy=    qenergy+vixen2*(1.0D0-vixen3/6.0D0)
C               IF (MGBWATER) THEN
C                  alpha_ij=SQRT(alpha(a)*alpha(b))
C                  vixen6= DUMMY2/(4*alpha_ij**2)
C                  fgb= SQRT(DUMMY2 + alpha_ij**2*exp(-vixen6))
C                  qenergy= qenergy- (1.0D0-vixen3/6.0D0)*0.5*(1-78**-1)*q(a)*q(b)/fgb
C               END IF

            ENDIF
         END DO
         DO c=1,qlistcount
            a=qlist(c,1)
            b=qlist(c,2)
            vixen4=q(a)
            DUMMY2=(x(a)-x(b))**2+(y(a)-y(b))**2+(z(a)-z(b))**2
            vixen5=SQRT(DUMMY2)
            NDUMMY1=1-bonds(a,b)
            NDUMMY3=NDUMMY1*(1-one_three(a,b))
            IF (NDUMMY3.NE.0) THEN
               vixen3=FLOAT(one_four(a,b))
               IF (FAKEWATER) THEN
                  vixen2=vixen4*q(b)/(dparam*dielec*vixen5**2) - vixen4*q(b)/(dparam*dielec*QCUTOFF**2)
               ELSE
                  vixen2=vixen4*q(b)/(vixen5*dielec*dparam) - vixen4*q(b)/(dparam*dielec*QCUTOFF)
               END IF
               qenergy=    qenergy+vixen2*(1.0D0-vixen3/6.0D0)
C               IF (MGBWATER) THEN
C                  alpha_ij=SQRT(alpha(a)*alpha(b))
C                  vixen6= DUMMY2/(4*alpha_ij**2)
C                  fgb= SQRT(DUMMY2 + alpha_ij**2*exp(-vixen6))
C                  qenergy= qenergy- (1.0D0-vixen3/6.0D0)*0.5*(1-78**-1)*q(a)*q(b)/fgb
C               END IF
            ENDIF
         END DO
      ELSE
         DO a=1,atoms
            vixen4=q(a)
            DO b=a+1,atoms
               DUMMY2=(x(a)-x(b))**2+(y(a)-y(b))**2+(z(a)-z(b))**2
C               IF (DUMMY2 .GT. AMCUTOFF**2) GOTO 31
               vixen5=SQRT(DUMMY2)
               r(a,b)=vixen5
               r(b,a)=vixen5

               NDUMMY1=1-bonds(a,b)
               NDUMMY3=NDUMMY1*(1-one_three(a,b))
               IF (NDUMMY3.NE.0) THEN
                  vixen3=FLOAT(one_four(a,b))

                  DUMMY=1.0D0/DUMMY2**3
                  vixen1=(vdwa(a,b)*DUMMY-vdwb(a,b))*DUMMY 
                  IF (FAKEWATER) THEN
                     vixen2=vixen4*q(b)/(dparam*dielec*vixen5**2) 
                  ELSE
                     vixen2=vixen4*q(b)/(vixen5*dielec*dparam) 
                  END IF

                  vdwenergy=vdwenergy+vixen1*(1.0D0-vixen3/2.0D0)
                  qenergy=    qenergy+vixen2*(1.0D0-vixen3/6.0D0)
                  IF (MGBWATER) THEN
                     alpha_ij=SQRT(alpha(a)*alpha(b))
                     vixen6= DUMMY2/(4*alpha_ij**2)
                     fgb= SQRT(DUMMY2 + alpha_ij**2*exp(-vixen6))
                     qenergy= qenergy- (1.0D0-vixen3/6.0D0)*0.5*(1-78**(-1))*q(a)*q(b)/fgb
                  END IF
               ENDIF
C31             CONTINUE
            ENDDO
         ENDDO
      END IF

      RETURN
      END


      SUBROUTINE chiraltest(CTEST,XVEC)
      USE commons
      USE modamber
      IMPLICIT NONE

      INTEGER  chiratom(4)
      INTEGER  marvin, J1
      DOUBLE PRECISION  tempxh,tempyh,tempzh,tempxi,tempyi,tempzi,tempxj,tempyj,tempzj
      DOUBLE PRECISION  tempxk,tempyk,tempzk,tempxl,tempyl,tempzl
      DOUBLE PRECISION  newxj,newyj,newzj,newxk,newyk,newzk,agamma,delta,dot
      DOUBLE PRECISION  cross(3)
      DOUBLE PRECISION  XVEC(3*NATOMS)
      LOGICAL  CTEST
      LOGICAL  move(NATOMS)
      INTEGER canine
      COMMON /CHIR/ canine

      CTEST=.FALSE.

      DO J1=1,atoms
         x(J1)=XVEC(3*J1-2)
         y(J1)=XVEC(3*J1-1)
         z(J1)=XVEC(3*J1)
      END DO
      DO a=1,atoms
        IF (chiral(a).NE.1) GOTO 17
        marvin=1
        DO b=1,atoms
          IF (bonds(a,b).EQ.1 .OR. bonds(b,a).EQ.1) THEN
            chiratom(marvin)=b
            marvin=marvin+1
          END IF
        END DO    

        i=chiratom(1)
        j=chiratom(2)
        k=chiratom(3)
        l=chiratom(4)

        IF (typech(i).EQ.'HC' .OR. typech(j).EQ.'HC' .OR. typech(k).EQ.'HC' .OR. typech(l).EQ.'HC') THEN
C
C  We are dealing with the beta carbon of isoleucine
C
        ELSE IF (typech(i).EQ.'OH'.OR.typech(j).EQ.'OH'.OR.typech(k).EQ.'OH'.OR.typech(l).EQ.'OH') THEN
C
C  We are dealing with the beta carbon of threonine
C  Note that we order the atoms backwards as this carbon should be 'R' configuration
C
C          PRINT *,' Found threonine beta carbon'
          IF (typech(i).EQ.'OH') THEN
            chiratom(4)=i
          ELSE IF (typech(j).EQ.'OH') THEN
            chiratom(4)=j
          ELSE IF (typech(k).EQ.'OH') THEN
            chiratom(4)=k
          ELSE
            chiratom(4)=l
          END IF
          IF (typech(i).EQ.'H1') THEN
            chiratom(1)=i
          ELSE IF (typech(j).EQ.'H1') THEN
            chiratom(1)=j
          ELSE IF (typech(k).EQ.'H1') THEN
            chiratom(1)=k
          ELSE
            chiratom(1)=l
          END IF
          IF (label(i).EQ.'d') THEN
            chiratom(3)=i
          ELSE IF (label(j).EQ.'d') THEN
            chiratom(3)=j
          ELSE IF (label(k).EQ.'d') THEN
            chiratom(3)=k
          ELSE 
            chiratom(3)=l
          END IF
          IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.j
     1       .OR. chiratom(3).EQ.j .OR. chiratom(4).EQ.j) .AND. (chiratom(1).EQ.k .OR. chiratom(3).EQ.k
     2       .OR. chiratom(4).EQ.k)) THEN
            chiratom(2)=l
          ELSE IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.j
     1       .OR. chiratom(3).EQ.j .OR. chiratom(4).EQ.j) .AND. (chiratom(1).EQ.l .OR. chiratom(3).EQ.l
     2       .OR. chiratom(4).EQ.l)) THEN          
            chiratom(2)=k
          ELSE IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.k
     1       .OR. chiratom(3).EQ.k .OR. chiratom(4).EQ.k) .AND. (chiratom(1).EQ.l .OR. chiratom(3).EQ.l
     2       .OR. chiratom(4).EQ.l)) THEN          
            chiratom(2)=j
          ELSE 
            chiratom(2)=i
          END IF

C          PRINT *,'Ranking for threonine beta carbon'
C          PRINT *,'1',typech(chiratom(4)),' ',chiratom(4)
C          PRINT *,'2',typech(chiratom(3)),' ',chiratom(3)
C          PRINT *,'3',typech(chiratom(2)),' ',chiratom(2)
C          PRINT *,'4',typech(chiratom(1)),' ',chiratom(1)

        ELSE
C  We have an alpha carbon
          IF (typech(j)(1:1).EQ.'H') THEN
            chiratom(1)=j
          ELSE IF (typech(k)(1:1).EQ.'H') THEN
            chiratom(1)=k
          ELSE IF (typech(l)(1:1).EQ.'H') THEN
            chiratom(1)=l
          END IF
    
          IF ((typech(i).EQ.'N ').OR.(typech(i).EQ.'N3')) THEN
            chiratom(2)=i
          ELSE IF ((typech(k).EQ.'N ').OR.(typech(k).EQ.'N3')) THEN
            chiratom(2)=k
          ELSE IF ((typech(l).EQ.'N ').OR.(typech(l).EQ.'N3')) THEN
            chiratom(2)=l
          END IF
  
          IF (typech(i).EQ.'C ') THEN
            chiratom(3)=i
          ELSE IF (typech(j).EQ.'C ') THEN
            chiratom(3)=j
          ELSE IF (typech(l).EQ.'C ') THEN
            chiratom(3)=l
          END IF
  
          IF (typech(i).EQ.'CT') THEN
            chiratom(4)=i
          ELSE IF (typech(j).EQ.'CT') THEN
            chiratom(4)=j
          ELSE IF (typech(k).EQ.'CT') THEN
            chiratom(4)=k
          END IF
        END IF

          h=a
          i=chiratom(1)
          j=chiratom(2)
          k=chiratom(3)
          l=chiratom(4)
C
C Move group to origin
C
        tempxi=x(i)-x(h)
        tempyi=y(i)-y(h)
        tempzi=z(i)-z(h)
        tempxj=x(j)-x(h)
        tempyj=y(j)-y(h)
        tempzj=z(j)-z(h)
        tempxk=x(k)-x(h)
        tempyk=y(k)-y(h)
        tempzk=z(k)-z(h)
        tempxl=x(l)-x(h)
        tempyl=y(l)-y(h)
        tempzl=z(l)-z(h)
        tempxh=0.0
        tempyh=0.0
        tempzh=0.0
  
        agamma=-(tempxj*tempxi+tempyj*tempyi+tempzj*tempzi)/(tempxi**2+tempyi**2+tempzi**2)
        delta=-(tempxk*tempxi+tempyk*tempyi+tempzk*tempzi)/(tempxi**2+tempyi**2+tempzi**2)       

        newxj=tempxj+agamma*tempxi
        newyj=tempyj+agamma*tempyi
        newzj=tempzj+agamma*tempzi
        newxk=tempxk+delta*tempxi
        newyk=tempyk+delta*tempyi
        newzk=tempzk+delta*tempzi
  
        cross(1)=newyj*newzk-newzj*newyk
        cross(2)=newzj*newxk-newxj*newzk
        cross(3)=newxj*newyk-newyj*newxk
  
        dot=cross(1)*tempxi+cross(2)*tempyi+cross(3)*tempzi
        IF (dot .LT. 0.0D0) THEN
C          PRINT *,'Configuration is S'
        ELSE
          CTEST=.TRUE. 
          WRITE(*,*) 'Configuration at atom ',a,' is wrong - dot= ',dot
          OPEN (UNIT=19, FILE='chiral.xyz')
          WRITE (19,*) natoms
          WRITE (19,*) ' '
          DO i=1,natoms
             WRITE (19,'(A1,3F20.10)') typech(i)(1:1),x(i),y(i),z(i)
          END DO
          CLOSE(19)
C
C Fix chiral centre - swap positions of side chain and H atom
C
          tempxi=x(chiratom(4)) - x(chiratom(1))
          tempyi=y(chiratom(4)) - y(chiratom(1))
          tempzi=z(chiratom(4)) - z(chiratom(1))
          x(chiratom(1))= x(chiratom(1))+tempxi
          y(chiratom(1))= y(chiratom(1))+tempyi
          z(chiratom(1))= z(chiratom(1))+tempzi
          DO J1=chiratom(4),atoms
             IF (bondedto(J1).GT.J1) GOTO 23
             IF (J1.NE.chiratom(1)) THEN
                WRITE (*,'(A,I4,A,3F8.3,A,3F8.3)') 'Moving atom ',J1,' from ',x(J1),y(J1),z(J1),' to ',x(J1)-tempxi,
     1          y(J1)-tempyi,z(J1)-tempzi
                x(J1)=x(J1)-tempxi
                y(J1)=y(J1)-tempyi
                z(J1)=z(J1)-tempzi
             END IF
          END DO 
23        CONTINUE
          OPEN (UNIT=19, FILE='fixed.xyz')
          WRITE (19,*) natoms
          WRITE (19,*) ' '
          DO i=1,natoms
             WRITE (19,'(A1,3F20.10)') typech(i)(1:1),x(i),y(i),z(i)
          END DO
          CLOSE(19)
        END IF
17    CONTINUE
      END DO

      DO J1=1,atoms
         XVEC(3*J1-2)=x(J1)
         XVEC(3*J1-1)=y(J1)
         XVEC(3*J1)=z(J1)
      END DO
  
      RETURN
      END

      SUBROUTINE TWIST(atom1,atom2,twistangle,NP)
      USE commons
      USE modamber
      IMPLICIT NONE
      
      INTEGER           atom1,atom2,atom3,NP
      DOUBLE PRECISION  twistangle,sa,ca,sp,cp,st,ct,m11,m12,m13,m21,m22,m23,m31,m32,m33,temptheta,tempphi
      DOUBLE PRECISION  s2t,c2t,s2p,c2p,oldx,oldy,oldz,newx,newy,newz


      DO a=1,atoms
        x(a)=COORDS(3*a-2,NP)
        y(a)=COORDS(3*a-1,NP)
        z(a)=COORDS(3*a,NP)
      END DO

C
C  Calculate elements of rotation matrix
C
      r(atom1,atom2)=SQRT((x(atom1)-x(atom2))**2+(y(atom1)-y(atom2))**2+(z(atom1)-z(atom2))**2)
C      PRINT *,'x=',x(atom2)-x(atom1)
C      PRINT *,'y=',y(atom2)-y(atom1)
C      PRINT *,'z=',z(atom2)-z(atom1)
C      PRINT *,' '
      temptheta=ACOS((z(atom2)-z(atom1))/r(atom1,atom2))
C      PRINT *,'Theta= ',temptheta*57.29577951
      st=SIN(temptheta)
      ct=COS(temptheta)
      s2t=st**2
      c2t=ct**2

      tempphi=ACOS((x(atom2)-x(atom1))/(r(atom1,atom2)*st))
      IF ((y(atom2)-y(atom1)).LT.0) tempphi=-tempphi
C      PRINT *,'Phi= ',tempphi*57.29577951
      sp=SIN(tempphi)
      cp=COS(tempphi)
      s2p=sp**2
      c2p=cp**2

      sa=SIN(twistangle)
      ca=COS(twistangle)
C
C  Rotation matrix elements    m11  m12  m13
C                              m21  m22  m23
C                              m31  m32  m33
C
      m11=(ca*c2t*c2p)+(ca*s2p)+(s2t*c2p)
      m12=(ca*c2t*sp*cp)+(sa*ct)-(ca*sp*cp)+(s2t*sp*cp)
      m13=(-ca*st*ct*cp)-(sa*st*sp)+(st*ct*cp)

      m21=(ca*c2t*sp*cp)-(sa*ct)-(ca*sp*cp)+(s2t*sp*cp)
      m22=(ca*c2t*s2p)+(ca*c2p)+(s2t*s2p)
      m23=(-ca*st*ct*sp)+(sa*st*cp)+(st*ct*sp)

      m31=(-ca*st*ct*cp)+(sa*st*sp)+(st*ct*cp)
      m32=(-ca*st*ct*sp)-(sa*st*cp)+(st*ct*sp)
      m33=(ca*s2t)+(c2t)
C
C  Decide whether it is a psi or a phi or an omega we are changing
C
      IF (typech(atom1).EQ.'N ') THEN
C
C  Phi
C
        DO a=atom2+1,atoms
          oldx=x(a)
          oldy=y(a)
          oldz=z(a)
C
C  Move system to origin (atom1)
C
          oldx=oldx-x(atom1)
          oldy=oldy-y(atom1)
          oldz=oldz-z(atom1)
C
C  Rotate
C
          newx = oldx*m11 + oldy*m12 + oldz*m13
          newy = oldx*m21 + oldy*m22 + oldz*m23
          newz = oldx*m31 + oldy*m32 + oldz*m33
C
C  Move back
C
          newx=newx+x(atom1)
          newy=newy+y(atom1)
          newz=newz+z(atom1)

          x(a)=newx
          y(a)=newy
          z(a)=newz
        END DO
      ELSE IF (typech(atom1).EQ.'CT') THEN
C
C  Psi
C
C  We need to find the next N
C
        DO a=atom2+1,atoms
          IF (typech(a).EQ.'N ' .AND. bonds(atom2,a).EQ.1) THEN
            atom3=a
            GOTO 117
          END IF
          IF (a.EQ.atoms) THEN
            PRINT *,'Unable to find next N in chain after atom',atom2
            STOP
          END IF
        END DO
117     CONTINUE
C
C  Move O in peptide link
C
        b=atom2+1
        oldx=x(b)
        oldy=y(b)
        oldz=z(b)

        oldx=oldx-x(atom1)
        oldy=oldy-y(atom1)
        oldz=oldz-z(atom1)

        newx = oldx*m11 + oldy*m12 + oldz*m13
        newy = oldx*m21 + oldy*m22 + oldz*m23
        newz = oldx*m31 + oldy*m32 + oldz*m33

        newx=newx+x(atom1)
        newy=newy+y(atom1)
        newz=newz+z(atom1)

        x(b)=newx
        y(b)=newy
        z(b)=newz
C
C  Now move the rest
C
        DO a=atom3-1,atoms
          oldx=x(a)
          oldy=y(a)
          oldz=z(a)
C
C  Move system to origin (atom1)
C
          oldx=oldx-x(atom1)
          oldy=oldy-y(atom1)
          oldz=oldz-z(atom1)
C
C  Rotate
C
          newx = oldx*m11 + oldy*m12 + oldz*m13
          newy = oldx*m21 + oldy*m22 + oldz*m23
          newz = oldx*m31 + oldy*m32 + oldz*m33
C
C  Move back
C
          newx=newx+x(atom1)
          newy=newy+y(atom1)
          newz=newz+z(atom1)

          x(a)=newx
          y(a)=newy
          z(a)=newz
        END DO

      ELSE IF (typech(atom1).EQ.'C ') THEN
CC
CC  Omega
CC
CC We need to find the next N
CC
C        DO b=atom1,atoms
C          IF (typech(b).NE.'N ') GOTO 209
C          atom3=b
C          GOTO 210
C209       CONTINUE
C        END DO

C210     CONTINUE
C        DO b=atom3-1,atoms
C          oldx=x(a)
C          oldy=y(a)
C          oldz=z(a)
CC
CC  Move system to origin (atom1)
CC
C          oldx=oldx-x(atom1)
C          oldy=oldy-y(atom1)
C          oldz=oldz-z(atom1)
CC
CC  Rotate
CC
C          newx = oldx*m11 + oldy*m12 + oldz*m13
C          newy = oldx*m21 + oldy*m22 + oldz*m23
C          newz = oldx*m31 + oldy*m32 + oldz*m33
CC
CC  Move back
CC
C          newx=newx+x(atom1)
C          newy=newy+y(atom1)
C          newz=newz+z(atom1)
C
C          x(a)=newx
C          y(a)=newy
C          z(a)=newz
C        END DO
C
      ELSE
        WRITE (*,*) 'Incorrect dihedral chosen :  ',typech(atom1),' - ',typech(atom2)
        STOP
      END IF

      DO a=1,atoms
        COORDS(3*a-2,NP)=x(a)
        COORDS(3*a-1,NP)=y(a)
        COORDS(3*a,NP)=z(a)
      END DO

      RETURN
      END
C
C  Step routine for amber.
C
      SUBROUTINE TAKESTEPAM(NP)
      USE commons
      USE modamber
      IMPLICIT NONE

      DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      COMMON /AMBWORD/    PMAX,PMIN,NMAX,NMIN,SIDESTEP

      DOUBLE PRECISION    P,ANGLE,DPRAND,RANDOM

      INTEGER             NP,NTEST1,NTEST2
C     CHARACTER  label(NATOMS)

C     CHARACTER(LEN=2) typech(NATOMS)

C     LOGICAL WATERSTEP
C     COMMON /WATER1/ WATERSTEP

      DO a=1,3*NATOMS
        COORDSO(a,NP)=COORDS(a,NP)
      ENDDO

C      OPEN (UNIT=12,FILE='prestep.xyz',STATUS='UNKNOWN')
C      WRITE (UNIT=12,FMT='(I4)') NATOMS
C      WRITE (UNIT=12,FMT='(A)') ' '
C      DO a=1,3*NATOMS
C        WRITE (UNIT=12,FMT='(A,3X,F14.10,3X,F14.10,3X,F14.10)') typech(a)(1:1),COORDS(3*a-2,NP),COORDS(3*a-1,NP),COORDS(3*a,NP)
C      END DO
C      CLOSE (UNIT=12,STATUS='KEEP')

      IF (WATERSTEP) THEN
         DO a=1,natoms
            IF (typech(a).EQ.'OW' .OR. typech(a).EQ.'HW') THEN
               RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
               COORDS(3*a-2,NP)=COORDS(3*a-2,NP)+RANDOM
               RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
               COORDS(3*a-1,NP)=COORDS(3*a-1,NP)+RANDOM
               RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
               COORDS(3*a,NP)=COORDS(3*a,NP)+RANDOM
            END IF
         END DO
      ELSE
192      CONTINUE
         b=0
         DO a=1,NDIHEDRALS
C
C  Calculate P, the probability of twisting
C
C           PRINT *,'a=',a
C           PRINT *,'NDIHEDRALS=',NDIHEDRALS
           IF (REAL(a).LE.(0.5*NDIHEDRALS)) THEN
             P=PMAX-a*((PMAX-PMIN)/(NDIHEDRALS*0.5))
           ELSE
             P=PMIN+(a-0.5*NDIHEDRALS)*((PMAX-PMIN)/(NDIHEDRALS*0.5))
           END IF
C           PRINT *,'P=',P

           ANGLE=(DPRAND()-0.5)*2.0*STEP(NP)
C           PRINT *,'ANGLE=',ANGLE
           RANDOM=DPRAND()
C           PRINT *,'RANDOM=',RANDOM
           IF (RANDOM.LT.P) THEN
              WRITE (*,'(A,I3,A,F10.5)') 'Twisting dihedral ',a,' by ',ANGLE*57.29577951
             CALL TWIST(DATOM1(a),DATOM2(a),ANGLE,NP)
             b=b+1
           END IF

         END DO

C        PRINT*,'NMIN,NMAX,NDIHEDRALS=',NMIN,NMAX,NDIHEDRALS
         NTEST1=INT(NMIN*NDIHEDRALS)
         IF (NTEST1.LT.1) NTEST1=1
         NTEST2=INT(NMAX*NDIHEDRALS)

C         WRITE (*,'(A,I3,A,I3,A)') 'Must shift between ',NTEST1,' and ',NTEST2,' dihedrals'
C         WRITE (*,'(A,I3)') 'Attempting to shift ',b

         IF (b.LT.NTEST1 .OR. b.GT.NTEST2) THEN
            WRITE (*,'(A)') 'Too many dihedrals shifted - retrying'
           GOTO 192
         END IF

C         OPEN (UNIT=12,FILE='midstep.xyz',STATUS='UNKNOWN')
C         WRITE (UNIT=12,FMT='(I4)') NATOMS
C         WRITE (UNIT=12,FMT='(A)') ' '
C         DO a=1,3*NATOMS
C           WRITE (UNIT=12,FMT='(A,3X,F14.10,3X,F14.10,3X,F14.10)') typech(a)(1:1),COORDS(3*a-2,NP),COORDS(3*a-1,NP),COORDS(3*a,NP)
C         END DO
C         CLOSE (UNIT=12,STATUS='KEEP')

          WRITE (*,'(A,I3,A,I3)') 'Shifting ',b,' dihedrals out of ',NDIHEDRALS
C
C Now shift sidechains by a small, random, cartesian displacement
C
         DO a=1,NATOMS
           IF (label(a).EQ.'d') GOTO 193
           RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
           COORDS(3*a-2,NP)=COORDS(3*a-2,NP)+RANDOM
           RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
           COORDS(3*a-1,NP)=COORDS(3*a-1,NP)+RANDOM
           RANDOM=(DPRAND()-0.5)*2.0*SIDESTEP
           COORDS(3*a,NP)=COORDS(3*a,NP)+RANDOM
193      CONTINUE
         END DO

C         OPEN (UNIT=12,FILE='poststep.xyz',STATUS='UNKNOWN')
C         WRITE (UNIT=12,FMT='(I4)') NATOMS
C         WRITE (UNIT=12,FMT='(A)') ' '
C         DO a=1,3*NATOMS
C           WRITE (UNIT=12,FMT='(A,3X,F14.10,3X,F14.10,3X,F14.10)') typech(a)(1:1),COORDS(3*a-2,NP),COORDS(3*a-1,NP),COORDS(3*a,NP)
C         END DO
C         CLOSE (UNIT=12,STATUS='KEEP')
      END IF


      RETURN

      END


      SUBROUTINE straighten(NP)
      USE commons
      USE modamber
      IMPLICIT NONE

      INTEGER            NP
      DOUBLE PRECISION   xcross,ycross,zcross,alph,temptheta,tempphi,ca,sa,ct,st,cp,sp,c2t,s2t,c2p,s2p
      DOUBLE PRECISION   oldx,oldy,oldz,m11,m12,m13,m21,m22,m23,m31,m32,m33,dot

      DO a=1,NATOMS
        x(a)=COORDS(3*a-2,NP)
        y(a)=COORDS(3*a-1,NP)
        z(a)=COORDS(3*a,NP)
      END DO

      OPEN (UNIT=7,FILE='before.xyz',STATUS='UNKNOWN')
      WRITE (UNIT=7,FMT='(I3)') atoms
      WRITE (UNIT=7,FMT='(A)') ' '
      DO b=1,atoms
        WRITE (UNIT=7,FMT='(A1,2X,F14.10,3X,F14.10,3X,F14.10)') typech(b)(1:1),x(b),y(b),z(b)
      END DO
      CLOSE (UNIT=7,STATUS='KEEP')

      DO a=1,t
        i=da1(a)
        j=da2(a)
        k=da3(a)
        l=da4(a)
        IF (label(j).NE.'d' .OR. label(k).NE.'d') GOTO 307
        IF (typech(i).NE.'O ' .OR. typech(j).NE.'C ' .OR. typech(k).NE.'N ' .OR. typech(l).NE.'H') GOTO 307

C        PRINT *,'Peptide link - atoms:'
C        WRITE (*,FMT='(I3,1X,I3,1X,I3,1X,I3)') i,j,k,l

        colin=0
        ax=x(j)-x(i)
        ay=y(j)-y(i)
        az=z(j)-z(i)
        bx=x(k)-x(j)
        by=y(k)-y(j)
        bz=z(k)-z(j)
        cx=x(j)-x(l)
        cy=y(j)-y(l)
        cz=z(j)-z(l)
        dx=x(l)-x(k)
        dy=y(l)-y(k)
        dz=z(l)-z(k)

        IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) .AND. ay*bz.GT.(by*az-TINY)
     1       .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
        IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) .AND. by*dz.GT.(dy*bz-TINY)
     1       .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2

        IF (colin.EQ.1) THEN
          PRINT *,'Three sites',i,j,k,'are colinear'
        ELSE IF (colin.EQ.2) THEN
          PRINT *,'Three sites',j,k,l,'are colinear'
        ELSE
          lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
          mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
          qx(1)=x(i)+(lambda*bx)
          qy(1)=y(i)+(lambda*by)
          qz(1)=z(i)+(lambda*bz)
          qx(2)=x(j)
          qy(2)=y(j)
          qz(2)=z(j)
          qx(3)=x(l)+(mu*bx)
          qy(3)=y(l)+(mu*by)
          qz(3)=z(l)+(mu*bz)
          numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
          denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1        *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))

          IF ((numer/denom).LT. -1) numer=numer+1D-12
          IF ((numer/denom).GT.1) numer=numer-1D-12

          dphi(a)=ACOS(numer/denom)
        ENDIF

        xcross=(qy(1)-qy(2))*(qz(3)-qz(2))-(qz(1)-qz(2))*(qy(3)-qy(2))
        ycross=(qz(1)-qz(2))*(qx(3)-qx(2))-(qx(1)-qx(2))*(qz(3)-qz(2))
        zcross=(qx(1)-qx(2))*(qy(3)-qy(2))-(qy(1)-qy(2))*(qx(3)-qx(2))

        dot=xcross*bx+ycross*by+zcross*bz
        IF (dot.LT.0.0D0) THEN
          alph=3.141592654-dphi(a)
        ELSE
          alph=dphi(a)-3.141592654
        END IF

C        PRINT *,'Rotating all atoms past',l,' by',alph*57.29577951
        ca=COS(alph)
        sa=SIN(alph)

        temptheta=ACOS(bz/SQRT(bx**2+by**2+bz**2))
C        PRINT *,'theta=',temptheta*57.29577951
        tempphi=ACOS(bx/(SQRT(bx**2+by**2+bz**2)*SIN(temptheta)))
        IF (by.LT.0.0D0) tempphi=-tempphi
C        PRINT *,'phi=',tempphi*57.29577951
        ct=COS(temptheta)
        st=SIN(temptheta)
        cp=COS(tempphi)
        sp=SIN(tempphi)

        c2t=ct**2
        c2p=cp**2
        s2t=st**2
        s2p=sp**2

C
C  Rotation matrix elements    m11  m12  m13
C                              m21  m22  m23
C                              m31  m32  m33
C
        m11=(ca*c2t*c2p)+(ca*s2p)+(s2t*c2p)
        m12=(ca*c2t*sp*cp)+(sa*ct)-(ca*sp*cp)+(s2t*sp*cp)
        m13=(-ca*st*ct*cp)-(sa*st*sp)+(st*ct*cp)

        m21=(ca*c2t*sp*cp)-(sa*ct)-(ca*sp*cp)+(s2t*sp*cp)
        m22=(ca*c2t*s2p)+(ca*c2p)+(s2t*s2p)
        m23=(-ca*st*ct*sp)+(sa*st*cp)+(st*ct*sp)

        m31=(-ca*st*ct*cp)+(sa*st*sp)+(st*ct*cp)
        m32=(-ca*st*ct*sp)-(sa*st*cp)+(st*ct*sp)
        m33=(ca*s2t)+(c2t)

        DO b=l,atoms
          x(b)=x(b)-x(j)
          y(b)=y(b)-y(j)
          z(b)=z(b)-z(j)

          oldx=x(b)
          oldy=y(b)
          oldz=z(b)

          x(b) = m11*oldx + m12*oldy + m13*oldz
          y(b) = m21*oldx + m22*oldy + m23*oldz
          z(b) = m31*oldx + m32*oldy + m33*oldz

          x(b)=x(b)+x(j)
          y(b)=y(b)+y(j)
          z(b)=z(b)+z(j)
        END DO

        DO b=1,atoms
          COORDS(3*b-2,NP)=x(b)
          COORDS(3*b-1,NP)=y(b)
          COORDS(3*b,NP)=z(b)
        END DO

307     CONTINUE
      END DO

      OPEN (UNIT=8,FILE='after.xyz',STATUS='UNKNOWN')
      WRITE (UNIT=8,FMT='(I3)') atoms
      WRITE (UNIT=8,FMT='(A)') ' '
      DO b=1,atoms
        WRITE (UNIT=8,FMT='(A1,2X,F14.10,3X,F14.10,3X,F14.10)') typech(b)(1:1),x(b),y(b),z(b)
      END DO
      CLOSE (UNIT=8,STATUS='KEEP')

      RETURN
      END


      SUBROUTINE ambermass
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER m1

      DO m1=1,atoms
         IF (typech(m1).EQ.'C ') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CA') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CB') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CC') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CK') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CM') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CN') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CQ') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CR') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CT') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CV') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CW') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'C*') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'C0') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'F ') THEN
            mass(m1)= 19.00
         ELSE IF (typech(m1).EQ.'H ') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HC') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H1') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H2') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H3') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HA') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H4') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H5') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HO') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HS') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HW') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HP') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'N ') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NA') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NB') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NC') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N2') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N3') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N*') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'O ') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OW') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OH') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OS') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'O2') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'P ') THEN
            mass(m1)= 30.97
         ELSE IF (typech(m1).EQ.'S ') THEN
            mass(m1)= 32.06
         ELSE IF (typech(m1).EQ.'SH') THEN
            mass(m1)= 32.06
         END IF
      END DO


      RETURN

      END 


      SUBROUTINE read_seedfile
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER  J1, nres

      WRITE (*,*) 'AMBERSEED= ',AMBERSEED
      OPEN (UNIT=27, FILE='amber.seed', STATUS='OLD')
      READ (27,*) nres
      DO J1=1,nres
         READ (27,*) seedphi(J1), seedpsi(J1)
         WRITE (*,'(A,I3,A,2F20.10)') 'Seed angles for residue ',J1,' are: ',seedphi(J1), seedpsi(J1)
         seedphi(J1)=seedphi(J1)/57.29577951
         seedpsi(J1)=seedpsi(J1)/57.29577951
      END DO
      CLOSE (27)

      RETURN

      END



      SUBROUTINE get_dihedrals
      USE commons
      USE modamber
      IMPLICIT NONE

      DOUBLE PRECISION torsvec1(3),torsvec2(3),crossproduct(3)
      DOUBLE PRECISION dotproduct
      INTEGER dummycount,dummycount2,dummycount3

      dummycount=1
      dummycount2=1
      dummycount3=1
      DO a=1,t
        i=da1(a)
        j=da2(a)
        k=da3(a)
        l=da4(a)
  
        IF ((label(i).NE.'d' .AND. label(i).NE.'t') .OR. label(j).NE.'d' .OR. label(k).NE.'d' .OR. 
     1     (label(l).NE.'d' .AND. label(l).NE.'t')) GOTO 13

        ax=x(j)-x(i)
        ay=y(j)-y(i)
        az=z(j)-z(i)
        bx=x(k)-x(j)
        by=y(k)-y(j)
        bz=z(k)-z(j)
        cx=x(j)-x(l)
        cy=y(j)-y(l)
        cz=z(j)-z(l)
        dx=x(l)-x(k)
        dy=y(l)-y(k)
        dz=z(l)-z(k)

        lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
        mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
        qx(1)=x(i)+(lambda*bx)
        qy(1)=y(i)+(lambda*by)
        qz(1)=z(i)+(lambda*bz)
        qx(2)=x(j)
        qy(2)=y(j)
        qz(2)=z(j)
        qx(3)=x(l)+(mu*bx)
        qy(3)=y(l)+(mu*by)
        qz(3)=z(l)+(mu*bz)
  
        numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
        denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1         *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))
   
        IF ((numer/denom).LT. -1) numer=numer+1D-12
        IF ((numer/denom).GT.1) numer=numer-1D-12
  
        dphi(a)=ACOS(numer/denom)

        torsvec1(1)=qx(1)-qx(2)
        torsvec1(2)=qy(1)-qy(2)
        torsvec1(3)=qz(1)-qz(2)
        torsvec2(1)=qx(3)-qx(2)
        torsvec2(2)=qy(3)-qy(2)
        torsvec2(3)=qz(3)-qz(2)
  
        crossproduct(1)=torsvec1(2)*torsvec2(3)-torsvec2(2)*torsvec1(3)
        crossproduct(2)=torsvec2(1)*torsvec1(3)-torsvec1(1)*torsvec2(3)
        crossproduct(3)=torsvec1(1)*torsvec2(2)-torsvec2(1)*torsvec1(2)

        actualphi(1)=0.0D0
        IF (dphi(a).LT.TINY .OR. dphi(a).GT.(3.14159)) THEN
          PRINT *,'Too close to 0 or 180 degrees - angle= ',dphi(a)
          IF ((typech(i).EQ.'N ' .OR. typech(i).EQ.'N3') .AND. typech(l).EQ.'N ') THEN
            IF (dummycount.EQ.1) dummycount2=2
            actualpsi(dummycount)=dphi(a)
            psiatom1(dummycount)=j
            psiatom2(dummycount)=k
            dummycount=dummycount+1
          ELSE IF (typech(i).EQ.'C ') THEN
            actualphi(dummycount2)=dphi(a)
            phiatom1(dummycount2)=j
            phiatom2(dummycount2)=k
            dummycount2=dummycount2+1
          ELSE IF (typech(i).EQ.'CT') THEN
            actualomega(dummycount3)=dphi(a)
            dummycount3=dummycount3+1
          ELSE
            PRINT *,'Wrong label - typech=',typech(i)
            STOP
          END IF
        ELSE
          dotproduct=bx*crossproduct(1)+by*crossproduct(2)+bz*crossproduct(3)
          IF (dotproduct.LT.0) THEN
            IF ((typech(i).EQ.'N ' .OR. typech(i).EQ.'N3') .AND. typech(l).EQ.'N ') THEN
              IF (dummycount.EQ.1) dummycount2=2
              psiatom1(dummycount)=j
              psiatom2(dummycount)=k
              actualpsi(dummycount)=-dphi(a)
              dummycount=dummycount+1
            ELSE IF (typech(i).EQ.'C ') THEN
              phiatom1(dummycount2)=j
              phiatom2(dummycount2)=k
              actualphi(dummycount2)=-dphi(a)
              dummycount2=dummycount2+1
            ELSE IF (typech(i).EQ.'CT') THEN
              actualomega(dummycount3)=-dphi(a)
              dummycount3=dummycount3+1
            ELSE
              PRINT *,'Wrong label - typech=',typech(i)
              STOP
            END IF          
          ELSE
            IF ((typech(i).EQ.'N ' .OR. typech(i).EQ.'N3') .AND. typech(l).EQ.'N ') THEN
              IF (dummycount.EQ.1) dummycount2=2
              actualpsi(dummycount)=dphi(a)
              psiatom1(dummycount)=j
              psiatom2(dummycount)=k
              dummycount=dummycount+1
            ELSE IF (typech(i).EQ.'C ') THEN
              phiatom1(dummycount2)=j
              phiatom2(dummycount2)=k
              actualphi(dummycount2)=dphi(a)
              dummycount2=dummycount2+1
            ELSE IF (typech(i).EQ.'CT') THEN
              actualomega(dummycount3)=dphi(a)
              dummycount3=dummycount3+1
            ELSE
              PRINT *,'Wrong label - typech=',typech(i)
              STOP
            END IF
         END IF    
        END IF
13    CONTINUE
      END DO

      res=dummycount-1
      DO i=1,res
         WRITE (*,'(A,I3,A,2F20.10)') 'Actual dihedral angles for residue ',i,' are: ',actualphi(i)*57.29577951, 
     1         57.29577951*actualpsi(i)
      END DO
      i=res+1
      IF (.NOT.(CAP)) WRITE (*,'(A,I3,A,2F20.10)') 'Actual dihedral angles for residue ',i,' are: ',actualphi(i)*57.29577951,
     1         57.29577951*actualpsi(i)
      
      RETURN
      END


      SUBROUTINE amseed
      USE commons
      USE modamber
      IMPLICIT NONE

      DOUBLE PRECISION  tempangle

      DO i=1,natoms
         COORDS(3*i-2,1)=x(i)
         COORDS(3*i-1,1)=y(i)
         COORDS(3*i,1)=z(i)
      END DO

      IF (CAP) THEN
         tempangle= actualphi(1)-seedphi(1)
         CALL twist(phiatom1(1),phiatom2(1),tempangle,1)
      END IF

      tempangle= actualpsi(1)-seedpsi(1)
      CALL twist(psiatom1(1),psiatom2(1),tempangle,1)

      DO i=2,res
         tempangle= actualphi(i)-seedphi(i)
         CALL twist(phiatom1(i),phiatom2(i),tempangle,1)
         tempangle= actualpsi(i)-seedpsi(i)
         CALL twist(psiatom1(i),psiatom2(i),tempangle,1)
      END DO

      IF (.NOT.(CAP)) THEN
         tempangle= actualphi(res+1)-seedphi(res+1)
         CALL twist(phiatom1(res+1),phiatom2(res+1),tempangle,1)
      END IF 

      OPEN (UNIT=27, FILE='seed.xyz')
      WRITE (27,*) natoms
      WRITE (27,*) ' '
      DO i=1,natoms
         WRITE (27,'(A,3F20.10)') typech(i)(1:1),x(i),y(i),z(i)
      END DO
      CLOSE(27)
       
      RETURN

      END

      SUBROUTINE GETALPHA
      USE commons
      USE modamber
      IMPLICIT NONE

      DO a=1,atoms
         IF (typech(a).EQ.'C ') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CA') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CB') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CC') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CK') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CM') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CN') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CQ') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CR') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CT') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CV') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'CW') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'C*') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'C0') THEN
            alpha(a)= 0.9615
         ELSE IF (typech(a).EQ.'F ') THEN
            alpha(a)= 1.0000
         ELSE IF (typech(a).EQ.'H ') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HC') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'H1') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'H2') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'H3') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HA') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'H4') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'H5') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HO') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HS') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HW') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'HP') THEN
            alpha(a)= 0.8461
         ELSE IF (typech(a).EQ.'N ') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'NA') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'NB') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'NC') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'N2') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'N3') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'N*') THEN
            alpha(a)= 0.9343
         ELSE IF (typech(a).EQ.'O ') THEN
            alpha(a)= 1.0088
         ELSE IF (typech(a).EQ.'OW') THEN
            alpha(a)= 1.0088
         ELSE IF (typech(a).EQ.'OH') THEN
            alpha(a)= 1.0088
         ELSE IF (typech(a).EQ.'OS') THEN
            alpha(a)= 1.0088
         ELSE IF (typech(a).EQ.'O2') THEN
            alpha(a)= 1.0088
         ELSE IF (typech(a).EQ.'P ') THEN
            alpha(a)= 1.0700
         ELSE IF (typech(a).EQ.'S ') THEN
            alpha(a)= 1.1733
         ELSE IF (typech(a).EQ.'SH') THEN
            alpha(a)= 1.1733
         END IF
      END DO

      RETURN

      END


      SUBROUTINE makelist
      USE commons
      USE modamber
      IMPLICIT NONE


      INTEGER  J1, J2
      DOUBLE PRECISION  DP1

C      WRITE (*,'(A)') 'Updating cutoff list'
      qlistcount=0
      rlistcount=0
      DO J1=1,atoms
         DO J2=J1+1,atoms
            IF (bonds(J1,J2).EQ.1 .OR. one_three(J1,J2).EQ.1) GOTO 139
            DP1=DSQRT((x(J1)-x(J2))**2 + (y(J1)-y(J2))**2 + (z(J1)-z(J2))**2)
            IF (DP1.LT.RCUTOFF) THEN
C               WRITE (*,*) 'Atom1, Atom2, R= ',J1,J2,DP1
               rlistcount= rlistcount+1
               rlist(rlistcount,1)=J1
               rlist(rlistcount,2)=J2
            ELSE IF (DP1.GE.RCUTOFF .AND. DP1.LT.QCUTOFF) THEN
C               WRITE (*,*) 'Atom1, Atom2, R= ',J1,J2,DP1
               qlistcount= qlistcount+1
               qlist(qlistcount,1)=J1
               qlist(qlistcount,2)=J2
            END IF
139         CONTINUE
         END DO
      END DO
C      WRITE (*,*) 'qlistcount= ',qlistcount
C      WRITE (*,*) 'rlistcount= ',rlistcount

      RETURN

      END

      SUBROUTINE moviedump
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER  J1
      LOGICAL  CRAP
      CHARACTER(LEN=1) adummy

      CRAP=.FALSE.
      WRITE (adummy,'(I1)') NQ(1)
      filename='movie.dump.' // adummy
      INQUIRE (FILE=filename, EXIST=CRAP)
      IF (CRAP) THEN
C        OPEN (UNIT=21,FILE=filename,ACCESS='APPEND') ! not allowed by Sun f90
         OPEN (UNIT=21,FILE=filename)
      ELSE
         OPEN (UNIT=21,FILE=filename)
      END IF
      WRITE (21,'(I6)') atoms
      WRITE (21,*) ' '
      DO J1=1,atoms
         WRITE (21,'(A1,3F20.10)') typech(J1)(1:1), x(J1), y(J1), z(J1)
      END DO
      CLOSE (21)

      RETURN
      END
      SUBROUTINE amberenergy
      USE commons
      USE modamber 
      IMPLICIT NONE
      DOUBLE PRECISION vixen1, vixen2, vixen3, vixen4, vixen5, vixen6, e1, e2, e3
      INTEGER canine
      COMMON /CHIR/ canine

C 
C First initialising count variables
C 
      benergy=0.0D0
      tenergy=0.0D0
      penergy=0.0D0
      vdwenergy=0.0D0
      totenergy=0.0D0
      impenergy=0.0D0
      qenergy=0.0D0
C 
C Calculating r,benergy,vdwenergy,qenergy - note that we calculate ALL r as we need them for vdW energy
C 
      CALL AMBERD

      DO i=1,bondnumber
         a=bondarray(i,1)
         b=bondarray(i,2)
         vixen1=DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
         benergy=benergy+kr(type(a),type(b))*((vixen1-ro(type(a),type(b)))**2)
C         WRITE (*,*) 'Bond ',a,' - ',b,'  kr= ',kr(type(a),type(b))
C         WRITE (*,*) 'Bond ',a,' - ',b,'  energy= ',kr(type(a),type(b))*((vixen1-ro(type(a),type(b)))**2)
      ENDDO

      CALL AMBERA

      DO a=1,ang
         b=aa1(a)
         c=aa2(a)
         d=aa3(a)

         vixen1=kt(type(b),type(c),type(d))*((theta(a)-to(type(b),type(c),type(d)))**2)
C         IF (chiral(c).EQ.1 .AND. ARMST) vixen1= vixen1*1000.0D0
         tenergy=tenergy+vixen1
C        PRINT *,'Angle ',b,c,d,typech(b),'-',typech(c),'-',typech(d),' theta=',theta(a)*57.29577951,' energy=',vixen1
      ENDDO

      DO a=1,t
         i=da1(a)
         j=da2(a)
         k=da3(a)
         l=da4(a)
         vixen1=dvn(a)/did(a)
         vixen2=dn(a)*dphi(a)-ddelta(a)
 
         vixen3=dvn2(a)/did(a)
         vixen4=dn2(a)*dphi(a)-ddelta2(a)

         vixen5=dvn3(a)/did(a)
         vixen6=dn3(a)*dphi(a)-ddelta3(a)

         penergy=penergy+(vixen1*(1.0D0+cos(vixen2)))+(vixen3*(1.0D0+cos(vixen4)))+(vixen5*(1.0D0+cos(vixen6)))
         E1=vixen1*(1.0D0+cos(vixen2))
         E2=(vixen3*(1.0D0+cos(vixen4)))
         E3=(vixen5*(1.0D0+cos(vixen6)))

C         PRINT *,'Torsion',i,j,k,l,typech(i),'-',typech(j),'-',typech(k),'-',typech(l),' Phi=',dphi(a)*57.29577951,
C     1           ' Energy=',E1+E2+E3
C         PRINT *,'PK=',dvn(a),' IDIVF=',did(a),' PN=',dn(a),' PHASE=',ddelta(a)

      ENDDO

      DO a=1,imp
         impenergy=impenergy+(ivn(a)*(1.0D0+COS(in1(a)*iphi(a)-idelta(a))))
     
C        PRINT *,'Improper',e,f,g,h,' energy= ',ivn(a)*(1+COS(in1(a)*iphi(a)-idelta(a)))

      ENDDO

      vixen2=0.0D0
C      DO a=1,atoms
C         IF (typech(a).EQ.'OW' .OR. typech(a).EQ.'HW') THEN
C            vixen1= DSQRT((x(a)-xbar)**2 + (y(a)-ybar)**2 + (z(a)-zbar)**2)
C         END IF
C         IF (vixen1.GT.30.0D0) THEN
C            WRITE (*,*) a,' vixen1= ',vixen1
C            vixen2= vixen2+ 1.0D4*(vixen1-30.0D0)**2
C         END IF
C      END DO

      IF (count.EQ.1) THEN 
         PRINT *,'ang= ',ang
         PRINT *,'t= ',t
         PRINT *,'imp= ',imp
      ENDIF
C     DO a=1,imp
C        PRINT *,ia1(a),ia2(a),ia3(a),ia4(a)
C        PRINT *,'phi= ',impropers(a)%phi*57.29577951
C     ENDDO

C      WRITE (*,*) 'Bond lengths'
C       DO a=1,atoms
C        DO b=a,atoms
C         IF (bonds(a,b).EQ.1) WRITE (*,*) a,'-',b,'    ',r(a,b),'  ',typech(a),' ',typech(b)
C        ENDDO
C       ENDDO

C       PRINT *,'Bond angles'
C         DO a=1,ang
C          PRINT *,angles(a)%a1,'-',angles(a)%a2,'-',angles(a)%a3,' = ',(angles(a)%theta)*57.29577951
C         ENDDO

C       PRINT *,'Torsion angles'
C       DO a=1,t
C        PRINT *,da1(a),torsions(a)%a2,torsions(a)%a3,torsions(a)%a4,'    ',&
C                57.29577951*torsions(a)%phi
C       ENDDO
    
C      PRINT *,' '
C      WRITE (*,'(A,F20.12)') 'Total bond strain energy= ',benergy
C      WRITE (*,'(A,F20.12)') 'Total angle strain energy= ',tenergy
C      WRITE (*,'(A,F20.12)') 'Torsion angle energy= ',penergy
C      WRITE (*,'(A,F20.12)') 'Improper torsion energy= ',impenergy
C      WRITE (*,'(A,F20.12)') 'vdW energy= ',vdwenergy
C      WRITE (*,'(A,F20.12)') 'qenergy= ',qenergy
C      'boundary energy= ',vixen2
      totenergy=benergy+tenergy+penergy+vdwenergy+impenergy+qenergy+vixen2
C       WRITE (*,'(A,F20.12)') 'Total energy= ',totenergy
C      IF (count.EQ.1) PRINT *,'Energy calculated'

      RETURN
      END

      SUBROUTINE AMBG
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER J1, qpower
      DOUBLE PRECISION dpbydx,dpbydy,dpbydz,dubydx,dubydy,dubydz,dvbydx,dvbydy,dvbydz,u,v
      DOUBLE PRECISION vixen1, vixen2, vixen3,vixen4,vixen5,vixen6,vixen7,vixen8,qcoeff,qbondfactor
      DOUBLE PRECISION vixen9, alpha_ij, fgb, DP1, DP2, DP3
      INTEGER canine
      COMMON /CHIR/ canine
C 
C Initialise all derivatives
C 
      DO J1=1,atoms
         dbondEbydx(J1)=0.0D0
         dbondEbydy(J1)=0.0D0
         dbondEbydz(J1)=0.0D0
         dangEbydx(J1)=0.0D0
         dangEbydy(J1)=0.0D0
         dangEbydz(J1)=0.0D0
         dtorsEbydx(J1)=0.0D0
         dtorsEbydy(J1)=0.0D0
         dtorsEbydz(J1)=0.0D0
         dvdwEbydx(J1)=0.0D0
         dvdwEbydy(J1)=0.0D0
         dvdwEbydz(J1)=0.0D0
         dqEbydx(J1)=0.0D0
         dqEbydy(J1)=0.0D0
         dqEbydz(J1)=0.0D0
      ENDDO

C      PRINT *,'Co-ordinates are:'
C      WRITE (*,FMT='(A,3F11.6)') 'O',x(1),y(1),z(1)
C      WRITE (*,FMT='(A,3F11.6)') 'H',x(2),y(2),z(2)
C      WRITE (*,FMT='(A,3F11.6)') 'H',x(3),y(3),z(3)

       DO c=1,bondnumber
          a=bondarray(c,1)
          b=bondarray(c,2)
          vixen1= DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
  
          dbondEbydx(a)=dbondEbydx(a)+(2.0D0*kr(type(a),type(b))*(x(a)-x(b))*(1.0D0-(ro(type(a),type(b))/vixen1)))
          dbondEbydy(a)=dbondEbydy(a)+(2.0D0*kr(type(a),type(b))*(y(a)-y(b))*(1.0D0-(ro(type(a),type(b))/vixen1)))
          dbondEbydz(a)=dbondEbydz(a)+(2.0D0*kr(type(a),type(b))*(z(a)-z(b))*(1.0D0-(ro(type(a),type(b))/vixen1)))
          dbondEbydx(b)=dbondEbydx(b)+(2.0D0*kr(type(a),type(b))*(x(b)-x(a))*(1.0D0-(ro(type(a),type(b))/vixen1)))
          dbondEbydy(b)=dbondEbydy(b)+(2.0D0*kr(type(a),type(b))*(y(b)-y(a))*(1.0D0-(ro(type(a),type(b))/vixen1)))
          dbondEbydz(b)=dbondEbydz(b)+(2.0D0*kr(type(a),type(b))*(z(b)-z(a))*(1.0D0-(ro(type(a),type(b))/vixen1)))

       ENDDO

C  PRINT *,"Bond derivatives calculated"
C 
C First calculate energy gradient from central atom in angle
C 
      DO d=1,ang    
        a=aa1(d)
        b=aa2(d)
        c=aa3(d)
 
        vixen1=(x(b)**2+y(b)**2+z(b)**2)-(x(b)*x(c)+y(b)*y(c)+z(b)*z(c))
        u=((x(c)-x(b))*x(a)+(y(c)-y(b))*y(a)+(z(c)-z(b))*z(a))+vixen1
        vixen2=DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
        vixen3=DSQRT((x(c)-x(b))**2 + (y(c)-y(b))**2 + (z(c)-z(b))**2)
        v=vixen2*vixen3
 
       dubydx=2.0D0*x(b)-x(c)-x(a)
       dubydy=2.0D0*y(b)-y(c)-y(a)
       dubydz=2.0D0*z(b)-z(c)-z(a)
       dvbydx=((x(b)-x(a))*vixen3**2+(x(b)-x(c))*vixen2**2)/v
       dvbydy=((y(b)-y(a))*vixen3**2+(y(b)-y(c))*vixen2**2)/v
       dvbydz=((z(b)-z(a))*vixen3**2+(z(b)-z(c))*vixen2**2)/v
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       e=type(a)
       f=type(b)
       g=type(c)

C       IF (chiral(b).EQ.1 .AND. ARMST) THEN
C          dangEbydx(b)=dangEbydx(b)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
C          dangEbydy(b)=dangEbydy(b)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
C          dangEbydz(b)=dangEbydz(b)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       ELSE
          dangEbydx(b)=dangEbydx(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
          dangEbydy(b)=dangEbydy(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
          dangEbydz(b)=dangEbydz(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       END IF
C 
C Now from end atoms
C 
       dubydx=x(c)-x(b)
       dubydy=y(c)-y(b)
       dubydz=z(c)-z(b)
       dvbydx=(x(a)-x(b))*vixen3/vixen2
       dvbydy=(y(a)-y(b))*vixen3/vixen2
       dvbydz=(z(a)-z(b))*vixen3/vixen2
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

C       IF (chiral(b).EQ.1 .AND. ARMST) THEN
C          dangEbydx(a)=dangEbydx(a)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
C          dangEbydy(a)=dangEbydy(a)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
C          dangEbydz(a)=dangEbydz(a)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       ELSE
          dangEbydx(a)=dangEbydx(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
          dangEbydy(a)=dangEbydy(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
          dangEbydz(a)=dangEbydz(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       END IF


       dubydx=x(a)-x(b)
       dubydy=y(a)-y(b)
       dubydz=z(a)-z(b)
       dvbydx=(x(c)-x(b))*vixen2/vixen3
       dvbydy=(y(c)-y(b))*vixen2/vixen3
       dvbydz=(z(c)-z(b))*vixen2/vixen3
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

C       IF (chiral(b).EQ.1 .AND. ARMST) THEN
C          dangEbydx(c)=dangEbydx(c)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
C          dangEbydy(c)=dangEbydy(c)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
C          dangEbydz(c)=dangEbydz(c)+2.0D3*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       ELSE
          dangEbydx(c)=dangEbydx(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
          dangEbydy(c)=dangEbydy(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
          dangEbydz(c)=dangEbydz(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C       END IF

      ENDDO

C  PRINT *,"Angle derivatives calculated"

       IF (FAKEWATER) THEN
         qcoeff=2.0D0
         qpower=3
       ELSE
         qcoeff=1.0D0
         qpower=2
       END IF

       IF (AMCUT) THEN
          DO c=1,rlistcount
             a=rlist(c,1)
             b=rlist(c,2)
             DP1=DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
             DP2=DBLE(1.0D0-(one_four(a,b)/2.0D0))
             DP3=DBLE(1.0D0-(one_four(a,b)/6.0D0))
             vixen1=DP2*((-12.0D0*vdwa(a,b)/DP1**13)+(6.0D0*vdwb(a,b))/DP1**7)
             vixen2=qcoeff*q(a)*q(b)*DP3/(dparam*dielec*DP1**qpower)
             vixen3=(x(a)-x(b))/DP1
             vixen4=(y(a)-y(b))/DP1
             vixen5=(z(a)-z(b))/DP1
             vixen6=-vixen3
             vixen7=-vixen4
             vixen8=-vixen5
C             vixen6=(x(b)-x(a))/DP1
C             vixen7=(y(b)-y(a))/DP1
C             vixen8=(z(b)-z(a))/DP1
C             IF (one_four(a,b).EQ.1 .AND. one_three(a,b).NE.1) THEN
C               dvdwEbydx(a)=dvdwEbydx(a)+0.5D0*(vixen1*vixen3)
C               dvdwEbydy(a)=dvdwEbydy(a)+0.5D0*(vixen1*vixen4)
C               dvdwEbydz(a)=dvdwEbydz(a)+0.5D0*(vixen1*vixen5)
C               dvdwEbydx(b)=dvdwEbydx(b)+0.5D0*(vixen1*vixen6)
C               dvdwEbydy(b)=dvdwEbydy(b)+0.5D0*(vixen1*vixen7)
C               dvdwEbydz(b)=dvdwEbydz(b)+0.5D0*(vixen1*vixen8)
C               dqEbydx(a)=dqEbydx(a)-(1.0D0/1.2D0)*vixen2*vixen3
C               dqEbydy(a)=dqEbydy(a)-(1.0D0/1.2D0)*vixen2*vixen4
C               dqEbydz(a)=dqEbydz(a)-(1.0D0/1.2D0)*vixen2*vixen5
C               dqEbydx(b)=dqEbydx(b)-(1.0D0/1.2D0)*vixen2*vixen6
C               dqEbydy(b)=dqEbydy(b)-(1.0D0/1.2D0)*vixen2*vixen7
C               dqEbydz(b)=dqEbydz(b)-(1.0D0/1.2D0)*vixen2*vixen8
C             ELSE IF (one_three(a,b).EQ.1) THEN
C              DO NOTHING
C             ELSE
                dvdwEbydx(a)=dvdwEbydx(a)+(vixen1*vixen3)
                dvdwEbydy(a)=dvdwEbydy(a)+(vixen1*vixen4)
                dvdwEbydz(a)=dvdwEbydz(a)+(vixen1*vixen5)
                dvdwEbydx(b)=dvdwEbydx(b)+(vixen1*vixen6)
                dvdwEbydy(b)=dvdwEbydy(b)+(vixen1*vixen7)
                dvdwEbydz(b)=dvdwEbydz(b)+(vixen1*vixen8)
                dqEbydx(a)=dqEbydx(a)-(vixen2*vixen3)
                dqEbydy(a)=dqEbydy(a)-(vixen2*vixen4)
                dqEbydz(a)=dqEbydz(a)-(vixen2*vixen5)
                dqEbydx(b)=dqEbydx(b)-(vixen2*vixen6)
                dqEbydy(b)=dqEbydy(b)-(vixen2*vixen7)
                dqEbydz(b)=dqEbydz(b)-(vixen2*vixen8)
C             END IF
          END DO
          DO c=1,qlistcount
             a=qlist(c,1)
             b=qlist(c,2)
             DP1=DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
             DP3=DBLE(1.0D0-(one_four(a,b)/6.0D0))
             vixen2=DP3*qcoeff*q(a)*q(b)/(dparam*dielec*DP1**qpower)
             vixen3=(x(a)-x(b))/DP1
             vixen4=(y(a)-y(b))/DP1
             vixen5=(z(a)-z(b))/DP1
             vixen6=-vixen3
             vixen7=-vixen4
             vixen8=-vixen5
C             vixen6=(x(b)-x(a))/DP1
C             vixen7=(y(b)-y(a))/DP1
C             vixen8=(z(b)-z(a))/DP1
C             IF (one_four(a,b).EQ.1 .AND. one_three(a,b).NE.1) THEN
C               dqEbydx(a)=dqEbydx(a)-(1.0D0/1.2D0)*vixen2*vixen3
C               dqEbydy(a)=dqEbydy(a)-(1.0D0/1.2D0)*vixen2*vixen4
C               dqEbydz(a)=dqEbydz(a)-(1.0D0/1.2D0)*vixen2*vixen5
C               dqEbydx(b)=dqEbydx(b)-(1.0D0/1.2D0)*vixen2*vixen6
C               dqEbydy(b)=dqEbydy(b)-(1.0D0/1.2D0)*vixen2*vixen7
C               dqEbydz(b)=dqEbydz(b)-(1.0D0/1.2D0)*vixen2*vixen8
C               IF (MGBWATER) THEN
C                  alpha_ij= SQRT(alpha(a)*alpha(b))
C                  vixen9= DP1**2/(4*alpha_ij**2)
C                  fgb= SQRT(DP1**2 + alpha_ij**2*exp(-vixen9))
C                  dqEbydx(a)= dqEbydx(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen3*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                  dqEbydy(a)= dqEbydy(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen4*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                  dqEbydz(a)= dqEbydz(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen5*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                  dqEbydx(b)= dqEbydx(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen6*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                  dqEbydy(b)= dqEbydy(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen7*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                  dqEbydz(b)= dqEbydz(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen8*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C               END IF
C             ELSE IF (one_three(a,b).EQ.1) THEN
C              DO NOTHING
C             ELSE
                dqEbydx(a)=dqEbydx(a)-(vixen2*vixen3)
                dqEbydy(a)=dqEbydy(a)-(vixen2*vixen4)
                dqEbydz(a)=dqEbydz(a)-(vixen2*vixen5)
                dqEbydx(b)=dqEbydx(b)-(vixen2*vixen6)
                dqEbydy(b)=dqEbydy(b)-(vixen2*vixen7)
                dqEbydz(b)=dqEbydz(b)-(vixen2*vixen8)
C                IF (MGBWATER) THEN
C                   alpha_ij= SQRT(alpha(a)*alpha(b))
C                   vixen9= DP1**2/(4*alpha_ij**2)
C                   fgb= SQRT(DP1**2 + alpha_ij**2*exp(-vixen9))
C                   dqEbydx(a)= dqEbydx(a) + 0.25*q(a)*q(b)*vixen3*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                   dqEbydy(a)= dqEbydy(a) + 0.25*q(a)*q(b)*vixen4*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                   dqEbydz(a)= dqEbydz(a) + 0.25*q(a)*q(b)*vixen5*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                   dqEbydx(b)= dqEbydx(b) + 0.25*q(a)*q(b)*vixen6*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                   dqEbydy(b)= dqEbydy(b) + 0.25*q(a)*q(b)*vixen7*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                   dqEbydz(b)= dqEbydz(b) + 0.25*q(a)*q(b)*vixen8*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
C                END IF
C             END IF
          END DO
       ELSE
          DO a=1,atoms
            DO b=a+1,atoms
              IF (bonds(a,b).NE.1) THEN
                DP1=DSQRT((x(a)-x(b))**2 + (y(a)-y(b))**2 + (z(a)-z(b))**2)
C                IF (DP1 .GT. AMCUTOFF) GOTO 37
                vixen1=(-12.0D0*vdwa(a,b)/DP1**13)+(6.0D0*vdwb(a,b))/DP1**7
                vixen2=qcoeff*q(a)*q(b)/(dparam*dielec*DP1**qpower)
                vixen3=(x(a)-x(b))/DP1
                vixen4=(y(a)-y(b))/DP1
                vixen5=(z(a)-z(b))/DP1
                vixen6=(x(b)-x(a))/DP1
                vixen7=(y(b)-y(a))/DP1
                vixen8=(z(b)-z(a))/DP1

                IF (one_four(a,b).EQ.1 .AND. one_three(a,b).NE.1) THEN
                  dvdwEbydx(a)=dvdwEbydx(a)+0.5D0*(vixen1*vixen3)
                  dvdwEbydy(a)=dvdwEbydy(a)+0.5D0*(vixen1*vixen4)
                  dvdwEbydz(a)=dvdwEbydz(a)+0.5D0*(vixen1*vixen5)
                  dvdwEbydx(b)=dvdwEbydx(b)+0.5D0*(vixen1*vixen6)
                  dvdwEbydy(b)=dvdwEbydy(b)+0.5D0*(vixen1*vixen7)
                  dvdwEbydz(b)=dvdwEbydz(b)+0.5D0*(vixen1*vixen8)
                  dqEbydx(a)=dqEbydx(a)-(1.0D0/1.2D0)*vixen2*vixen3
                  dqEbydy(a)=dqEbydy(a)-(1.0D0/1.2D0)*vixen2*vixen4
                  dqEbydz(a)=dqEbydz(a)-(1.0D0/1.2D0)*vixen2*vixen5
                  dqEbydx(b)=dqEbydx(b)-(1.0D0/1.2D0)*vixen2*vixen6
                  dqEbydy(b)=dqEbydy(b)-(1.0D0/1.2D0)*vixen2*vixen7
                  dqEbydz(b)=dqEbydz(b)-(1.0D0/1.2D0)*vixen2*vixen8
                  IF (MGBWATER) THEN
                     alpha_ij= SQRT(alpha(a)*alpha(b))
                     vixen9= DP1**2/(4*alpha_ij**2)
                     fgb= SQRT(DP1**2 + alpha_ij**2*exp(-vixen9))
                     dqEbydx(a)= dqEbydx(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen3*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                     dqEbydy(a)= dqEbydy(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen4*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                     dqEbydz(a)= dqEbydz(a) + 0.25*q(a)*q(b)*(1/1.2)*vixen5*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                     dqEbydx(b)= dqEbydx(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen6*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                     dqEbydy(b)= dqEbydy(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen7*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                     dqEbydz(b)= dqEbydz(b) + 0.25*q(a)*q(b)*(1/1.2)*vixen8*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                  END IF

                ELSE IF (one_three(a,b).EQ.1) THEN
C                 DO NOTHING
                ELSE 
                   dvdwEbydx(a)=dvdwEbydx(a)+(vixen1*vixen3)
                   dvdwEbydy(a)=dvdwEbydy(a)+(vixen1*vixen4)
                   dvdwEbydz(a)=dvdwEbydz(a)+(vixen1*vixen5)
                   dvdwEbydx(b)=dvdwEbydx(b)+(vixen1*vixen6)
                   dvdwEbydy(b)=dvdwEbydy(b)+(vixen1*vixen7)
                   dvdwEbydz(b)=dvdwEbydz(b)+(vixen1*vixen8)
                   dqEbydx(a)=dqEbydx(a)-(vixen2*vixen3)
                   dqEbydy(a)=dqEbydy(a)-(vixen2*vixen4)
                   dqEbydz(a)=dqEbydz(a)-(vixen2*vixen5)
                   dqEbydx(b)=dqEbydx(b)-(vixen2*vixen6)
                   dqEbydy(b)=dqEbydy(b)-(vixen2*vixen7)
                   dqEbydz(b)=dqEbydz(b)-(vixen2*vixen8)
                   IF (MGBWATER) THEN
                      alpha_ij= SQRT(alpha(a)*alpha(b))
                      vixen9= DP1**2/(4*alpha_ij**2)
                      fgb= SQRT(DP1**2 + alpha_ij**2*exp(-vixen9))
                      dqEbydx(a)= dqEbydx(a) + 0.25*q(a)*q(b)*vixen3*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                      dqEbydy(a)= dqEbydy(a) + 0.25*q(a)*q(b)*vixen4*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                      dqEbydz(a)= dqEbydz(a) + 0.25*q(a)*q(b)*vixen5*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                      dqEbydx(b)= dqEbydx(b) + 0.25*q(a)*q(b)*vixen6*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                      dqEbydy(b)= dqEbydy(b) + 0.25*q(a)*q(b)*vixen7*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                      dqEbydz(b)= dqEbydz(b) + 0.25*q(a)*q(b)*vixen8*DP1*0.5*(1-78**(-1))*(4 - exp(-vixen9))/fgb**3
                   END IF

                 END IF
               END IF
37             CONTINUE
             ENDDO
           ENDDO
       END IF
211    CONTINUE

      CALL pderivs

      CALL impderivs
 
      DO a=1,atoms
       dEbydx(a)=dbondEbydx(a)+dangEbydx(a)+dvdwEbydx(a)+dtorsEbydx(a)+dimpEbydx(a)+dqEbydx(a)
       dEbydy(a)=dbondEbydy(a)+dangEbydy(a)+dvdwEbydy(a)+dtorsEbydy(a)+dimpEbydy(a)+dqEbydy(a)
       dEbydz(a)=dbondEbydz(a)+dangEbydz(a)+dvdwEbydz(a)+dtorsEbydz(a)+dimpEbydz(a)+dqEbydz(a)
      ENDDO

C      DO a=1,atoms
C         IF (typech(a).EQ.'OW' .OR. typech(a).eq.'HW') THEN
C            vixen1= DSQRT((x(a)-xbar)**2+(y(a)-ybar)**2+(z(a)-zbar)**2)
C            IF (vixen1.GT.30.0D0) THEN
C               dEbydx(a)= dEbydx(a)+ 2.0D4*(vixen1-30.0D0)*(x(a)-xbar)/vixen1
C               dEbydy(a)= dEbydy(a)+ 2.0D4*(vixen1-30.0D0)*(y(a)-ybar)/vixen1
C               dEbydz(a)= dEbydz(a)+ 2.0D4*(vixen1-30.0D0)*(z(a)-zbar)/vixen1
C            END IF
C         END IF
C      END DO

      IF (count.EQ.1) PRINT *,'Derivatives calculated'
      count=0

      RETURN
      END

      SUBROUTINE pderivs
      USE commons
      USE modamber
      IMPLICIT NONE
      DOUBLE PRECISION xe,xf,xg,ye,yf,yg,ze,zf,zg,u,v,up,vp,p
      DOUBLE PRECISION dlambdabydx,dlambdabydy,dlambdabydz,dmubydx,dmubydy,dmubydz
      DOUBLE PRECISION dpbydx,dpbydy,dpbydz,dubydx,dubydy,dubydz,dvbydx,dvbydy,dvbydz
      DOUBLE PRECISION dphibydx,dphibydy,dphibydz,dupbydx,dupbydy,dupbydz,dvpbydx,dvpbydy,dvpbydz
      DOUBLE PRECISION vixen1, vixen2, vixen3, dp1 
      INTEGER J1
C 
C Initialise derivatives
C 
      DO J1=1,atoms
         dtorsEbydx(J1)=0.0D0
         dtorsEbydy(J1)=0.0D0
         dtorsEbydz(J1)=0.0D0
      ENDDO
C 
C Calculate derivatives for all atoms in position "a"
C 
      DO j=1,t
C 
C Initialise intermediates
C 
       dubydx=0.0D0
       dubydy=0.0D0
       dubydz=0.0D0
       dvbydx=0.0D0
       dvbydy=0.0D0
       dvbydz=0.0D0
       dpbydx=0.0D0
       dpbydy=0.0D0
       dupbydx=0.0D0
       dupbydy=0.0D0
       dupbydz=0.0D0
       dvpbydx=0.0D0
       dvpbydy=0.0D0
       dvpbydz=0.0D0
       dpbydz=0.0D0

       i=da1(j)
       IF (i.EQ.0) GOTO 20
       a=da1(j)
       b=da2(j)
       c=da3(j)
       d=da4(j)
C 
C First define all functions as on paper
C 
       dp1=DSQRT((x(b)-x(c))**2 + (y(b)-y(c))**2 + (z(b)-z(c))**2)
       xe=x(c)-x(b)
       ye=y(c)-y(b)
       ze=z(c)-z(b)
       lambda=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(a)*xe)+(y(a)*ye)+(z(a)*ze)))/(xe**2+ye**2+ze**2)
       mu=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(d)*xe)+(y(d)*ye)+(z(d)*ze)))/(xe**2+ye**2+ze**2)
       xf=x(b)-x(a)-(lambda*xe)
       yf=y(b)-y(a)-(lambda*ye)
       zf=z(b)-z(a)-(lambda*ze)
       xg=x(b)-x(d)-(mu*xe)
       yg=y(b)-y(d)-(mu*ye)
       zg=z(b)-z(d)-(mu*ze)
       vixen1=(x(b)**2+y(b)**2+z(b)**2)+mu*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)-((x(b)*x(d))+(y(b)*y(d))+(z(b)*z(d)))
       vixen3=(x(d)-x(b))*x(a)+(y(d)-y(b))*y(a)+(z(d)-z(b))*z(a)
       u=vixen1+vixen2+vixen3+(lambda*mu*dp1**2)
       up=SQRT(xf**2+yf**2+zf**2)
       vp=SQRT(xg**2+yg**2+zg**2)
       v=up*vp
       p=u/v

       IF (p.LT.-1) p=p+TINY
       IF (p.GT.1) p=p-TINY
C 
C Now calculate derivatives
C 
       dlambdabydx=-xe/dp1**2
       dlambdabydy=-ye/dp1**2
       dlambdabydz=-ze/dp1**2
       dubydx=(mu*xe)+x(d)-x(b)+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))
       dubydy=(mu*ye)+y(d)-y(b)+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))
       dubydz=(mu*ze)+z(d)-z(b)+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))

       dvbydx=(vp/up)*(-xf-dlambdabydx*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydy=(vp/up)*(-yf-dlambdabydy*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydz=(vp/up)*(-zf-dlambdabydz*((xf*xe)+(yf*ye)+(zf*ze)))

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
C 
C Now add terms to energy gradients
C 
       CALL hairy(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position "b"
C 
       i=da2(j)      

       vixen1=(dp1**2)*(x(a)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(x(b)-x(c))
       dlambdabydx=vixen1/(dp1**4)
       vixen1=(dp1**2)*(y(a)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(y(b)-y(c))
       dlambdabydy=vixen1/(dp1**4)
       vixen1=(dp1**2)*(z(a)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(z(b)-z(c))
       dlambdabydz=vixen1/(dp1**4)

       vixen1=(dp1**2)*(x(d)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(b)-x(c))
       dmubydx=vixen1/(dp1**4)
       vixen1=(dp1**2)*(y(d)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(b)-y(c))
       dmubydy=vixen1/(dp1**4)
       vixen1=(dp1**2)*(z(d)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(b)-z(c))
       dmubydz=vixen1/(dp1**4)

       vixen1=mu*(2.0D0*x(b)-x(c)-x(a))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*x(b)-x(c)-x(d))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydx)+(mu*dlambdabydx))*dp1**2
       dubydx=2.0D0*x(b)-x(a)-x(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(x(b)-x(c))
  
       vixen1=mu*(2.0D0*y(b)-y(c)-y(a))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*y(b)-y(c)-y(d))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydy)+(mu*dlambdabydy))*dp1**2
       dubydy=2.0D0*y(b)-y(a)-y(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(y(b)-y(c))

       vixen1=mu*(2.0D0*z(b)-z(c)-z(a))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*z(b)-z(c)-z(d))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydz)+(mu*dlambdabydz))*dp1**2
       dubydz=2.0D0*z(b)-z(a)-z(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(z(b)-z(c))
       
       dupbydx=(xf*(lambda+1)-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(yf*(lambda+1)-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(zf*(lambda+1)-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(xg*(mu+1)-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(yg*(mu+1)-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(zg*(mu+1)-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
      
       CALL hairy(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position 3
C 
       i=da3(j)

       dlambdabydx=(dp1**2*(x(b)-x(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (x(c)-x(b)))/dp1**4    
       dlambdabydy=(dp1**2*(y(b)-y(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (y(c)-y(b)))/dp1**4
       dlambdabydz=(dp1**2*(z(b)-z(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (z(c)-z(b)))/dp1**4     

       dmubydx=(dp1**2*(x(b)-x(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(c)-x(b)))
     1         /dp1**4
       dmubydy=(dp1**2*(y(b)-y(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(c)-y(b)))
     1         /dp1**4
       dmubydz=(dp1**2*(z(b)-z(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(c)-z(b)))
     1         /dp1**4

       vixen1=mu*(x(a)-x(b))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydx)+(mu*dlambdabydx))*dp1**2+2.0D0*lambda*mu*xe
       dubydx=vixen1+lambda*(x(d)-x(b))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(y(a)-y(b))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydy)+(mu*dlambdabydy))*dp1**2+2.0D0*lambda*mu*ye
       dubydy=vixen1+lambda*(y(d)-y(b))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(z(a)-z(b))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydz)+(mu*dlambdabydz))*dp1**2+2.0D0*lambda*mu*ze
       dubydz=vixen1+lambda*(z(d)-z(b))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       dupbydx=(-lambda*xf-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(-lambda*yf-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(-lambda*zf-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(-mu*xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(-mu*yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(-mu*zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairy(dpbydx,dpbydy,dpbydz)       
C 
C Now for all atoms in position 4
C 
       i=da4(j)

       dmubydx=-xe/dp1**2
       dmubydy=-ye/dp1**2
       dmubydz=-ze/dp1**2

       vixen1=lambda*dmubydx*dp1**2
       vixen2=lambda*dmubydy*dp1**2
       vixen3=lambda*dmubydz*dp1**2
       dubydx=(lambda*xe)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydx)+x(a)-x(b)+vixen1
       dubydy=(lambda*ye)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydy)+y(a)-y(b)+vixen2
       dubydz=(lambda*ze)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydz)+z(a)-z(b)+vixen3

       dvbydx=up*(-xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvbydy=up*(-yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvbydz=up*(-zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairy(dpbydx,dpbydy,dpbydz)
      ENDDO 
20    CONTINUE

      RETURN
      END

      SUBROUTINE hairy(dpbydx,dpbydy,dpbydz)
      USE commons
      USE modamber
      IMPLICIT NONE
      DOUBLE PRECISION dpbydx,dpbydy,dpbydz
      DOUBLE PRECISION vixen1,vixen2,vixen3

       e=type(a)
       f=type(b)
       g=type(c)
       h=type(d)
       PK=dvn(j)
       PK2=dvn2(j)
       PK3=dvn3(j)
       PN=dn(j)
       PN2=dn2(j)
       PN3=dn3(j)
       PHASE=ddelta(j)
       PHASE2=ddelta2(j)
       PHASE3=ddelta3(j)
       IDIVF=did(j)

       vixen1=PK*PN/IDIVF
       vixen2=PK2*PN2/IDIVF
       vixen3=PK3*PN3/IDIVF       

       IF (PN.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE)
       ELSE IF (PN.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
       ELSE IF (PN.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
       ENDIF

       IF (PN2.EQ.1) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*COS(PHASE2)
       ELSE IF (PN2.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE2)
       ELSE IF (PN2.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
       ELSE IF (PN2.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
       ENDIF

       IF (PN3.EQ.1) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*COS(PHASE3)
       ELSE IF (PN3.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE3)
       ELSE IF (PN3.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
       ELSE IF (PN3.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
       ENDIF

       RETURN
       END

      SUBROUTINE hairyimp(dpbydx,dpbydy,dpbydz)
      USE commons
      USE modamber
      IMPLICIT NONE
      DOUBLE PRECISION dpbydx,dpbydy,dpbydz,vixen1

      e=type(a)
      f=type(b)
      g=type(c)
      h=type(d)
      IPK=ivn(j)
      IPHASE=idelta(j)
      IPN=in1(j)

      vixen1=IPK*IPN
   
      IF (IPN.EQ.2) THEN
        dimpEbydx(i)=dimpEbydx(i)+2.0D0*dpbydx*vixen1*COS(iphi(j))*COS(IPHASE)       
        dimpEbydy(i)=dimpEbydy(i)+2.0D0*dpbydy*vixen1*COS(iphi(j))*COS(IPHASE)       
        dimpEbydz(i)=dimpEbydz(i)+2.0D0*dpbydz*vixen1*COS(iphi(j))*COS(IPHASE)       
      ELSE IF (IPN.EQ.0) THEN

      ELSE 
        PRINT *,'IPN for this torsion is not catered for in derivatives'
        STOP
      ENDIF

      RETURN
      END
       
      SUBROUTINE impderivs
      USE commons
      USE modamber
      IMPLICIT NONE
      DOUBLE PRECISION xe,xf,xg,ye,yf,yg,ze,zf,zg,u,v,up,vp,p
      DOUBLE PRECISION dlambdabydx,dlambdabydy,dlambdabydz,dmubydx,dmubydy,dmubydz
      DOUBLE PRECISION dpbydx,dpbydy,dpbydz,dubydx,dubydy,dubydz,dvbydx,dvbydy,dvbydz
      DOUBLE PRECISION dphibydx,dphibydy,dphibydz,dupbydx,dupbydy,dupbydz,dvpbydx,dvpbydy,dvpbydz
      DOUBLE PRECISION vixen1, vixen2, vixen3, dp1
      INTEGER J1
C 
C Initialise derivatives
C 
      DO J1=1,atoms
         dimpEbydx(J1)=0.0D0
         dimpEbydy(J1)=0.0D0
         dimpEbydz(J1)=0.0D0
      ENDDO
C 
C Calculate derivatives for all atoms in position "a"
C 
      DO j=1,imp
C 
C Initialise intermediates
C 
       dubydx=0.0D0
       dubydy=0.0D0
       dubydz=0.0D0
       dvbydx=0.0D0
       dvbydy=0.0D0
       dvbydz=0.0D0
       dpbydx=0.0D0
       dpbydy=0.0D0
       dupbydx=0.0D0
       dupbydy=0.0D0
       dupbydz=0.0D0
       dvpbydx=0.0D0
       dvpbydy=0.0D0
       dvpbydz=0.0D0
       dpbydz=0.0D0

       i=ia1(j)
       IF (i.EQ.0) GOTO 10
       a=ia1(j)
       b=ia2(j)
       c=ia3(j)
       d=ia4(j)
C 
C First define all functions as on paper
C 
       dp1=DSQRT((x(c)-x(b))**2 + (y(c)-y(b))**2 + (z(c)-z(b))**2)
       xe=x(c)-x(b)
       ye=y(c)-y(b)
       ze=z(c)-z(b)
       lambda=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(a)*xe)+(y(a)*ye)+(z(a)*ze)))/(xe**2+ye**2+ze**2)
       mu=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(d)*xe)+(y(d)*ye)+(z(d)*ze)))/(xe**2+ye**2+ze**2)
       xf=x(b)-x(a)-(lambda*xe)
       yf=y(b)-y(a)-(lambda*ye)
       zf=z(b)-z(a)-(lambda*ze)
       xg=x(b)-x(d)-(mu*xe)
       yg=y(b)-y(d)-(mu*ye)
       zg=z(b)-z(d)-(mu*ze)
       vixen1=(x(b)**2+y(b)**2+z(b)**2)+mu*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)-((x(b)*x(d))+(y(b)*y(d))+(z(b)*z(d)))
       vixen3=(x(d)-x(b))*x(a)+(y(d)-y(b))*y(a)+(z(d)-z(b))*z(a)
       u=vixen1+vixen2+vixen3+(lambda*mu*dp1**2)
       up=SQRT(xf**2+yf**2+zf**2)
       vp=SQRT(xg**2+yg**2+zg**2)
       v=up*vp
       p=u/v

        IF (p.LT.-1) p=p+TINY
        IF (p.GT.1) p=p-TINY
C 
C Now calculate derivatives
C 
       dlambdabydx=-xe/dp1**2
       dlambdabydy=-ye/dp1**2
       dlambdabydz=-ze/dp1**2
       dubydx=(mu*xe)+x(d)-x(b)+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))
       dubydy=(mu*ye)+y(d)-y(b)+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))
       dubydz=(mu*ze)+z(d)-z(b)+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*dp1**2))

       dvbydx=(vp/up)*(-xf-dlambdabydx*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydy=(vp/up)*(-yf-dlambdabydy*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydz=(vp/up)*(-zf-dlambdabydz*((xf*xe)+(yf*ye)+(zf*ze)))

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
C 
C Now add terms to energy gradients
C 
       CALL hairyimp(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position "b"
C 
       i=ia2(j)      

       vixen1=(dp1**2)*(x(a)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(x(b)-x(c))
       dlambdabydx=vixen1/(dp1**4)
       vixen1=(dp1**2)*(y(a)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(y(b)-y(c))
       dlambdabydy=vixen1/(dp1**4)
       vixen1=(dp1**2)*(z(a)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(z(b)-z(c))
       dlambdabydz=vixen1/(dp1**4)

       vixen1=(dp1**2)*(x(d)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(b)-x(c))
       dmubydx=vixen1/(dp1**4)
       vixen1=(dp1**2)*(y(d)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(b)-y(c))
       dmubydy=vixen1/(dp1**4)
       vixen1=(dp1**2)*(z(d)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(b)-z(c))
       dmubydz=vixen1/(dp1**4)

       vixen1=mu*(2.0D0*x(b)-x(c)-x(a))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*x(b)-x(c)-x(d))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydx)+(mu*dlambdabydx))*dp1**2
       dubydx=2.0D0*x(b)-x(a)-x(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(x(b)-x(c))
  
       vixen1=mu*(2.0D0*y(b)-y(c)-y(a))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*y(b)-y(c)-y(d))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydy)+(mu*dlambdabydy))*dp1**2
       dubydy=2.0D0*y(b)-y(a)-y(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(y(b)-y(c))

       vixen1=mu*(2.0D0*z(b)-z(c)-z(a))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*z(b)-z(c)-z(d))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydz)+(mu*dlambdabydz))*dp1**2
       dubydz=2.0D0*z(b)-z(a)-z(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(z(b)-z(c))
       
       dupbydx=(xf*(lambda+1)-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(yf*(lambda+1)-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(zf*(lambda+1)-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(xg*(mu+1)-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(yg*(mu+1)-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(zg*(mu+1)-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
      
       CALL hairyimp(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position 3
C 
       i=ia3(j)

       dlambdabydx=(dp1**2*(x(b)-x(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (x(c)-x(b)))/dp1**4    
       dlambdabydy=(dp1**2*(y(b)-y(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (y(c)-y(b)))/dp1**4
       dlambdabydz=(dp1**2*(z(b)-z(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (z(c)-z(b)))/dp1**4     

       dmubydx=(dp1**2*(x(b)-x(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(c)-x(b)))
     1         /dp1**4
       dmubydy=(dp1**2*(y(b)-y(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(c)-y(b)))
     1         /dp1**4
       dmubydz=(dp1**2*(z(b)-z(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(c)-z(b)))
     1         /dp1**4

       vixen1=mu*(x(a)-x(b))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydx)+(mu*dlambdabydx))*dp1**2+2.0D0*lambda*mu*xe
       dubydx=vixen1+lambda*(x(d)-x(b))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(y(a)-y(b))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydy)+(mu*dlambdabydy))*dp1**2+2.0D0*lambda*mu*ye
       dubydy=vixen1+lambda*(y(d)-y(b))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(z(a)-z(b))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydz)+(mu*dlambdabydz))*dp1**2+2.0D0*lambda*mu*ze
       dubydz=vixen1+lambda*(z(d)-z(b))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       dupbydx=(-lambda*xf-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(-lambda*yf-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(-lambda*zf-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(-mu*xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(-mu*yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(-mu*zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairyimp(dpbydx,dpbydy,dpbydz)       
C 
C Now for all atoms in position 4
C 
       i=ia4(j)

       dmubydx=-xe/dp1**2
       dmubydy=-ye/dp1**2
       dmubydz=-ze/dp1**2

       vixen1=lambda*dmubydx*dp1**2
       vixen2=lambda*dmubydy*dp1**2
       vixen3=lambda*dmubydz*dp1**2
       dubydx=(lambda*xe)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydx)+x(a)-x(b)+vixen1
       dubydy=(lambda*ye)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydy)+y(a)-y(b)+vixen2
       dubydz=(lambda*ze)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydz)+z(a)-z(b)+vixen3

       dvbydx=up*(-xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvbydy=up*(-yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvbydz=up*(-zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairyimp(dpbydx,dpbydy,dpbydz)
      ENDDO 
10    CONTINUE

      RETURN
      END


      SUBROUTINE numergrad
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER i1
      DOUBLE PRECISION         numerderiv(3*NATOMS)
      DOUBLE PRECISION         et

      DO i1=1,3*NATOMS
         numerderiv(i1)= 0.0
      END DO

      DO i1=1,atoms
         x(i1)= x(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         x(i1)= x(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1-2)= (et-totenergy)/0.00002
         x(i1)= x(i1)+ 0.00001

         y(i1)= y(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         y(i1)= y(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1-1)= (et-totenergy)/0.00002
         y(i1)= y(i1)+ 0.00001

         z(i1)= z(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         z(i1)= z(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1)= (et-totenergy)/0.00002
         z(i1)= z(i1)+ 0.00001
      END DO

      WRITE (*,FMT='(8X,A,5X,A,5X,A,5X,A,5X,A)') 'NUMERICAL','   *   ','ANALYTIC','   *   ','DIFFERENCE'
      WRITE (*,FMT='(8X,A,5X,A,5X,A,5X,A,5X,A)') '---------','   *   ','--------','   *   ','----------'
      DO i1=1,atoms
         WRITE (*,FMT='(I3,F20.10,A,F20.10,A,F20.10)') i1,numerderiv(3*i1-2),'  * ',
     1         dEbydx(i1),'   * ',dEbydx(i1)-numerderiv(3*i1-2)
         WRITE (*,FMT='(I3,F20.10,A,F20.10,A,F20.10)') i1,numerderiv(3*i1-1),'  * ',
     1         dEbydy(i1),'   * ',dEbydy(i1)-numerderiv(3*i1-1)
         WRITE (*,FMT='(I3,F20.10,A,F20.10,A,F20.10)') i1,numerderiv(3*i1),'  * ',
     1         dEbydz(i1),'   * ',dEbydz(i1)-numerderiv(3*i1)
      END DO


      RETURN

      END

      SUBROUTINE aparams
      USE commons
      USE modamber
      IMPLICIT NONE
      DOUBLE PRECISION vixen1, vixen2, vixen3, vixen4
     
      PRINT *,'Reading parameters'

      OPEN (UNIT=9,FILE='amber.dat',STATUS='OLD')
C 
C Reading parameters
C 
C Format of params.dat for bonds  typechb typechc kr(typechb,typechc) ro(typechb,typechc)
C 
      DO a=1,200
         READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
         IF (ios.LT.0) GOTO 10
         IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') GOTO 10
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,1X,A2)') typechb,typechc

         CALL typenumber(typechb)
         b=ans
         CALL typenumber(typechc)
         c=ans 
         IF (b.EQ.-1 .OR. c.EQ.-1) THEN
          PRINT *,'Bond parameters ',typechb,'-',typechc,' specified incorrectly - aborting'
          STOP
         ENDIF

         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,1X,A2,2X,F5.0,4X,F6.0)') typechb,typechc,kr(b,c),ro(b,c)
C        WRITE(*,'(3I3,1X,A,2X,A,2X,2F15.5)') a,b,c,typechb,typechc,kr(b,c),ro(b,c)
         kr(c,b)=kr(b,c)
         ro(c,b)=ro(b,c)
      ENDDO
10    CONTINUE
C 
C Format of params.dat for angles  type1 type2 type3 kt(1,2,3) degto(1,2,3)
C 
      DO a=1,300
         READ (UNIT=9,IOSTAT=ios,FMT='(A5)') check 
         IF (ios.LT.0) GOTO 20
         IF (.NOT.(check.EQ.'start' .OR. check.EQ.'Start' .OR. check.EQ.'START')) THEN
         IF (.NOT.(check.NE.'start' .AND. check.NE.'Start' .AND. check.NE.'START' .AND. a.EQ.1)) THEN
         IF (check.EQ.'end' .OR. check.EQ.'End' .OR. check.EQ.'END') GOTO 20
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,1X,A2,1X,A2)') typechb,typechc,typechd

         CALL typenumber(typechb)
         b=ans
         CALL typenumber(typechc)
         c=ans
         CALL typenumber(typechd)
         d=ans

         IF (b.EQ.-1 .OR. c.EQ.-1 .OR. d.EQ.-1) THEN
            PRINT *,'Angle ',typechb,'-',typechc,'-',typechd,' specified incorrectly - aborting'
            STOP
         ENDIF
   
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,3X,F5.0,6X,F6.0)') typechb,typechc,typechd,kt(b,c,d),degto(b,c,d)
         to(b,c,d)=degto(b,c,d)/57.29577951
         kt(d,c,b)=kt(b,c,d)
         to(d,c,b)=to(b,c,d)
         degto(d,c,b)=degto(b,c,d)
         ENDIF
         ENDIF
      ENDDO
20    CONTINUE

      IF (a.EQ.1) THEN
       PRINT *,'Error reading angle parameters - no start'
       STOP
      ENDIF 
C 
C Format of params.dat for torsions  type1 type2 type3 type4 PK IDIVF PN PHASE
C 
      DO a=1,100
         READ (UNIT=9,IOSTAT=ios,FMT='(A5)') check
         IF (ios.LT.0) GOTO 30
         IF (.NOT.(check.EQ.'start' .OR. check.EQ.'Start' .OR. check.EQ.'START')) THEN
         IF (check.NE.'start' .AND. check.NE.'Start' .AND. check.NE.'START' .AND. a.EQ.1) GOTO 30
         IF (check.EQ.'end' .OR. check.EQ.'End' .OR. check.EQ.'END')  GOTO 30
         IF (a.GT.1) BACKSPACE 9
         
         READ (UNIT=9,IOSTAT=ios,FMT='(A2)') typechb
         IF (typechb.NE.'X ') GOTO 30
         BACKSPACE 9

         READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2,3X,F1.0,3X,F6.0,7X,F5.0,13X,F2.0)') typechb,typechc,typechd,
     1        typeche,vixen1,vixen2,vixen3,vixen4
         CALL typenumber(typechc)
           c=ans
         CALL typenumber(typechd)
           d=ans

         gentorsparams(a-1,1)=c
         gentorsparams(a-1,2)=d
         gentorsparams(a-1,3)=vixen1 
         gentorsparams(a-1,4)=vixen2
         gentorsparams(a-1,5)=vixen3/57.29577951
         gentorsparams(a-1,6)=ABS(vixen4)

         ENDIF
      ENDDO
30    CONTINUE
 
      IF (a.EQ.1) THEN
       PRINT *,'Error reading torsion parameters - no start'
      STOP
      ENDIF

      BACKSPACE 9
      DO a=1,50
        READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
C        PRINT *,"check= ",check
        IF (check.EQ."end") THEN
C        PRINT *,"End of torsion parameters reached - a=",a
        GOTO 32
        END IF
        BACKSPACE 9

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,1X,A2,1X,A2,1X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F3.0)') typechb,typechc,
     1       typechd,typeche,vixen1,vixen2,vixen3,vixen4
        CALL typenumber(typechb)
          b=ans
        CALL typenumber(typechc)
          c=ans
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans

        spectorsparams(a,1)=b
        spectorsparams(a,2)=c
        spectorsparams(a,3)=d
        spectorsparams(a,4)=e
        spectorsparams(a,5)=vixen1
        spectorsparams(a,6)=vixen2
        spectorsparams(a,7)=vixen3/57.29577951
        spectorsparams(a,8)=ABS(vixen4)

        IF (vixen4.GT.0) THEN
         spectorsparams(a,9)=0.0
         spectorsparams(a,10)=0.0
         spectorsparams(a,11)=0.0
         spectorsparams(a,12)=0.0
         spectorsparams(a,13)=0.0
         spectorsparams(a,14)=0.0
         GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,1X,A2,1X,A2,1X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') typechb,typechc,
     1      typechd,typeche,vixen1,vixen2,vixen3,vixen4

        spectorsparams(a,9)=vixen2
        spectorsparams(a,10)=vixen3/57.29577951
        spectorsparams(a,11)=ABS(vixen4)
       
        IF (vixen4.GT.0) THEN
          spectorsparams(a,12)=0.0
          spectorsparams(a,13)=0.0
          spectorsparams(a,14)=0.0
          GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,1X,A2,1X,A2,1X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') typechb,typechc,
     1    typechd,typeche,vixen1,vixen2,vixen3,vixen4

        spectorsparams(a,12)=vixen2
        spectorsparams(a,13)=vixen3/57.29577951
        spectorsparams(a,14)=vixen4

31     CONTINUE
       END DO
32     CONTINUE
C 
C Format for improper torsion parameters    typecha  typechb  typechc  typechd  IPK  IPHASE  IPN
C 
      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A9)') check
       IF (ios.LT.0) THEN
        PRINT *,'End of parameters file without improper torsion parameters'
        STOP
       ENDIF
       IF (.NOT.(check.NE.'impropers' .AND. check.NE.'IMPROPERS' .AND. check.NE.'Impropers')) GOTO 40
      ENDDO
40    CONTINUE

      IF (a.GE.100) THEN 
       PRINT *,'No improper torsion parameters specified'
      STOP
      ENDIF

      bean=0
      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A2,1X,A2,1X,A2,1X,A2)') typechb,typechc,typechd,typeche
       IF (ios.LT.0) THEN
        PRINT *,'End of file during improper torsion parameters'
      STOP
       ENDIF
      IF (typechb.EQ.'X ' .AND. typechc.EQ.'X ') THEN
        BACKSPACE 9 
        READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2)') typechb,typechc,typechd,typeche
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans
 
        genimpparams(a-bean,1)=d
        genimpparams(a-bean,2)=e

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3
        genimpparams(a-bean,3)=vixen1
        genimpparams(a-bean,4)=vixen2/57.29577951
        genimpparams(a-bean,5)=vixen3

       ELSE IF (typechb.EQ.'X ' .AND. typechc.NE.'X ') THEN
        bean=bean+1
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2)') typechb,typechc,typechd,typeche
         CALL typenumber(typechc)
           c=ans
         CALL typenumber(typechd)
           d=ans
         CALL typenumber(typeche)
           e=ans

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3

        midimpparams(bean,1)=c
        midimpparams(bean,2)=d        
        midimpparams(bean,3)=e
        midimpparams(bean,4)=vixen1
        midimpparams(bean,5)=vixen2/57.29577951
        midimpparams(bean,6)=vixen3

       ELSE IF (typechb.EQ.'  ') THEN
          GOTO 50

       ELSE IF (typechb.EQ.'en' .OR. typechb.EQ.'En' .OR. typechb.EQ.'EN') THEN
          GOTO 60
       ELSE 
          GOTO 60

       ENDIF
50     CONTINUE
      ENDDO   
60    CONTINUE

      BACKSPACE 9
      DO a=1,20
        READ (UNIT=9,FMT='(A2,1X,A2,1X,A2,1X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3
        IF (typechb.EQ.'en') GOTO 65
        CALL typenumber(typechb)
          b=ans
        CALL typenumber(typechc)
          c=ans
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans

        specimpparams(a,1)=b
        specimpparams(a,2)=c
        specimpparams(a,3)=d
        specimpparams(a,4)=e
        specimpparams(a,5)=vixen1
        specimpparams(a,6)=vixen2/57.29577951
        specimpparams(a,7)=vixen3

      END DO
65    CONTINUE

      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
       IF (ios.LT.0) THEN
        PRINT *,'End of parameter file without vdW parameters'
        STOP
       ENDIF
       IF (.NOT.(check.NE.'vdw' .AND. check.NE.'vdW')) GOTO 70
      ENDDO
70    CONTINUE
  
      IF (a.LT.100) THEN
       DO a=1,42
        READ (UNIT=9,IOSTAT=ios,FMT='(A2)') typechb
        IF (typechb.EQ.'en') GOTO 80
        IF (ios.LT.0) GOTO 80
        CALL typenumber(typechb)
        b=ans
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,10X,F6.4,2X,F6.4)') typechb,vdwr(b),vdwe(b)
       ENDDO
80     CONTINUE
      ELSE
       PRINT *,'Unable to find vdW parameters'
      ENDIF

      PRINT *,'Parameters read'
    
      RETURN
      END

      SUBROUTINE typenumber(typechar)
      USE modamber
      IMPLICIT NONE
      CHARACTER(LEN=2)  typechar

      IF (typechar.EQ.'C ') THEN
         ans=1
      ELSE IF (typechar.EQ.'CA') THEN
         ans=2
      ELSE IF (typechar.EQ.'CB') THEN
         ans=3
      ELSE IF (typechar.EQ.'CC') THEN
         ans=4
      ELSE IF (typechar.EQ.'CK') THEN
         ans=5
      ELSE IF (typechar.EQ.'CM') THEN
         ans=6
      ELSE IF (typechar.EQ.'CN') THEN
         ans=7
      ELSE IF (typechar.EQ.'CQ') THEN
         ans=8
      ELSE IF (typechar.EQ.'CR') THEN
         ans=9
      ELSE IF (typechar.EQ.'CT') THEN
         ans=10
      ELSE IF (typechar.EQ.'CV') THEN
         ans=11
      ELSE IF (typechar.EQ.'CW') THEN
         ans=12
      ELSE IF (typechar.EQ.'C*') THEN
         ans=13
      ELSE IF (typechar.EQ.'C0') THEN
         ans=14
      ELSE IF (typechar.EQ.'F ') THEN
         ans=15
      ELSE IF (typechar.EQ.'H ') THEN
         ans=16
      ELSE IF (typechar.EQ.'HC') THEN
         ans=17
      ELSE IF (typechar.EQ.'H1') THEN
         ans=18
      ELSE IF (typechar.EQ.'H2') THEN
         ans=19
      ELSE IF (typechar.EQ.'H3') THEN
         ans=20
      ELSE IF (typechar.EQ.'HA') THEN
         ans=21
      ELSE IF (typechar.EQ.'H4') THEN
         ans=22
      ELSE IF (typechar.EQ.'H5') THEN
         ans=23
      ELSE IF (typechar.EQ.'HO') THEN
         ans=24
      ELSE IF (typechar.EQ.'HS') THEN
         ans=25
      ELSE IF (typechar.EQ.'HW') THEN
         ans=26
      ELSE IF (typechar.EQ.'HP') THEN
         ans=27
      ELSE IF (typechar.EQ.'N ') THEN
         ans=28
      ELSE IF (typechar.EQ.'NA') THEN
         ans=29
      ELSE IF (typechar.EQ.'NB') THEN
         ans=30
      ELSE IF (typechar.EQ.'NC') THEN
         ans=31
      ELSE IF (typechar.EQ.'N2') THEN
         ans=32
      ELSE IF (typechar.EQ.'N3') THEN
         ans=33
      ELSE IF (typechar.EQ.'N*') THEN
         ans=34
      ELSE IF (typechar.EQ.'O ') THEN
         ans=35
      ELSE IF (typechar.EQ.'OW') THEN
         ans=36
      ELSE IF (typechar.EQ.'OH') THEN
         ans=37
      ELSE IF (typechar.EQ.'OS') THEN
         ans=38
      ELSE IF (typechar.EQ.'O2') THEN
         ans=39
      ELSE IF (typechar.EQ.'P ') THEN
         ans=40
      ELSE IF (typechar.EQ.'S ') THEN
         ans=41
      ELSE IF (typechar.EQ.'SH') THEN
         ans=42
      ELSE 
         ans=-1
      ENDIF

      RETURN
      END
      SUBROUTINE aread
      USE commons
      USE modamber
      IMPLICIT NONE
      INTEGER canine
      COMMON /CHIR/ canine
      DOUBLE PRECISION vixen1

      ALLOCATE(DATOM1(NATOMS),DATOM2(NATOMS))
      ALLOCATE(bondarray(NATOMS*4,2))
      ALLOCATE(x(NATOMS),y(NATOMS),z(NATOMS),q(NATOMS),mass(NATOMS))
      ALLOCATE(r(NATOMS,NATOMS),vdwa(NATOMS,NATOMS),vdwb(NATOMS,NATOMS))
      ALLOCATE(qlist(NATOMS**2,2),rlist(NATOMS**2,2))
      ALLOCATE(aa1(NATOMS*5),aa2(NATOMS*5),aa3(NATOMS*5))
      ALLOCATE(theta(NATOMS*5))
      ALLOCATE(da1(NATOMS*5),da2(NATOMS*5),da3(NATOMS*5),da4(NATOMS*5))
      ALLOCATE(dphi(NATOMS*5),did(NATOMS*5),dvn(NATOMS*5),dvn2(NATOMS*5))
      ALLOCATE(dvn3(NATOMS*5),ddelta(NATOMS*5),ddelta2(NATOMS*5),ddelta3(NATOMS*5))
      ALLOCATE(dn(NATOMS*5),dn2(NATOMS*5),dn3(NATOMS*5))
      ALLOCATE(ia1(NATOMS*5),ia2(NATOMS*5),ia3(NATOMS*5),ia4(NATOMS*5))
      ALLOCATE(iphi(NATOMS*5),ivn(NATOMS*5),idelta(NATOMS*5),in1(NATOMS*5))
      ALLOCATE(atnum(NATOMS),bondedto(NATOMS),type(NATOMS))
      ALLOCATE(bonds(NATOMS,NATOMS))
      ALLOCATE(one_four(NATOMS,NATOMS),one_three(NATOMS,NATOMS))
      ALLOCATE(label(NATOMS))
      ALLOCATE(typech(NATOMS))
      ALLOCATE(dbondEbydx(NATOMS),dbondEbydy(NATOMS),dbondEbydz(NATOMS),dangEbydx(NATOMS),
     &                 dangEbydy(NATOMS),dangEbydz(NATOMS),
     &                 dtorsEbydx(NATOMS),dtorsEbydy(NATOMS),dtorsEbydz(NATOMS),dvdwEbydx(NATOMS),
     &                 dvdwEbydy(NATOMS),dvdwEbydz(NATOMS),
     &                 dEbydx(NATOMS),dEbydy(NATOMS),dEbydz(NATOMS),dimpEbydx(NATOMS),
     &                 dimpEbydy(NATOMS),dimpEbydz(NATOMS),
     &                 dqEbydx(NATOMS),dqEbydy(NATOMS),dqEbydz(NATOMS))
C      ALLOCATE(bondhell(3*NATOMS,3*NATOMS), anglehell(3*NATOMS,3*NATOMS), torshell(3*NATOMS,3*NATOMS),
C     1         imphell(3*NATOMS,3*NATOMS), qhell(3*NATOMS,3*NATOMS), vdwhell(3*NATOMS,3*NATOMS),
C     2         hell(3*NATOMS,3*NATOMS))
      ALLOCATE(CHIRAL(NATOMS),CHIRALARRAY(NATOMS,5))

      OPEN (UNIT=9,FILE='coords.amber',STATUS='OLD')
C 
C Format for info file   label  typech  number_in_list  bondedto  x  y  z
C 
C      PRINT *,'Max. no. of atoms = ',NATOMS
      DO a=1,NATOMS
        READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
        IF (ios.LT.0) THEN
          PRINT *,'End of file before all information specified',a
          STOP
        END IF
        IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') GOTO 10
        BACKSPACE 9
        READ (UNIT=9,FMT=*) label(a),typech(a),atnum(a),bondedto(a),x(a),
     1       y(a),z(a)
        IF (label(a).EQ.'*') THEN
          chiral(a)=1
C         PRINT *,'Found chiral atom',a
        END IF
      END DO
10    CONTINUE

      PRINT *,'Co-ordinates read successfully'
      atoms=a-1
      PRINT *,'atoms= ',atoms
C 
C Check for any rings and close them!
C 
      DO a=1,10
       READ (UNIT=9,IOSTAT=ios,FMT=*) check
       IF (check.NE.'loop' .AND. a.EQ.10) GOTO 20
       IF (.NOT.(check.NE.'loop')) THEN
          BACKSPACE 9
          READ (UNIT=9,FMT='(A4,7X,I2)') check,rings
          GOTO 20
      ENDIF
      END DO
20    CONTINUE

      bean=1
      IF (a.LT.10) THEN
       DO i=1,rings
        READ (UNIT=9,FMT=*) loopatom(2*i-1),loopatom(2*i)
         b=loopatom(2*i-1)
         c=loopatom(2*i)
         IF (b.EQ.0 .OR. c.EQ.0) THEN
          PRINT *,'No atoms specified for loop bond'
         END IF
        bonds(b,c)=1
       bonds(c,b)=1
       bondarray(bean,1)=b
       bondarray(bean,2)=c
       bean=bean+1
       END DO
      ELSE
       PRINT *,'No loop line in co-ordinate file'
       STOP
      END IF
C 
C Read charges
C 
      DO a=1,10
       READ (UNIT=9,IOSTAT=ios,FMT=*) check
       IF (check.NE.'charges' .AND. a.EQ.10) GOTO 30
       IF (.NOT.(check.NE.'charges')) GOTO 30
      END DO
30    CONTINUE

      IF (a.EQ.10) THEN
       PRINT *,'No charges specified'
       STOP
      ELSE 
      DO c=1,atoms
       READ (UNIT=9,IOSTAT=ios,FMT=*) b,vixen1
C       WRITE (*,*) b,'  ',vixen1
       q(b)=vixen1*18.2223
       IF (ios.LT.0) GOTO 40
      END DO
40    CONTINUE
      END IF
C 
C Converting typech to type number
C 
      DO a=1,atoms
       CALL typenumber(typech(a))
       type(a)=ans
       IF (type(a).EQ.-1) THEN
        PRINT *,'Atom ',a,'type specified incorrectly - aborting'
      STOP
       END IF
      END DO
C 
C Connecting backbone atoms
C 
      DO a=2,atoms
       IF (label(a).NE.'b') GOTO 50
       bonds(a,a-1)=1
       bonds(a-1,a)=1
      END DO       
50    CONTINUE
C 
C Connecting side-chain atoms
C 
      DO a=2,atoms
       IF (.NOT.(label(a).EQ.'b')) THEN
        IF (bondedto(a).NE.0) THEN
         bonds(a,bondedto(a))=1
         bonds(bondedto(a),a)=1
         bondarray(bean,1)=a
         bondarray(bean,2)=bondedto(a)
         bean=bean+1
        END IF
       ENDIF
      END DO

      bondnumber=bean-1

      CALL one3

      RETURN
      END

      SUBROUTINE one3
      USE commons
      USE modamber
      IMPLICIT NONE
      PRINT *,'SUBROUTINE one3'
       DO a=1,atoms
         DO b=a+1,atoms
           IF (bonds(a,b) .NE. 1) THEN
             DO c=1,atoms
               IF (bonds(a,c).EQ.1 .AND. bonds(b,c).EQ.1) THEN
                 one_three(a,b)=1
                 one_three(b,a)=1
C                PRINT *,a,b
               END IF
             END DO
           END IF
         END DO
       END DO

      RETURN

      END
      SUBROUTINE AMBERA
      USE commons
      USE modamber
      IMPLICIT NONE

      DOUBLE PRECISION  vixen1

      DO a=1,ang
         b=aa1(a)
         c=aa2(a)
         d=aa3(a)

         vixen1=DSQRT((x(b)-x(c))**2 + (y(b)-y(c))**2 + (z(b)-z(c))**2)
         top=(x(b)-x(c))*(x(d)-x(c))+(y(b)-y(c))*(y(d)-y(c))+(z(b)-z(c))*(z(d)-z(c))
         bottom=vixen1*SQRT((x(d)-x(c))**2+(y(d)-y(c))**2+(z(d)-z(c))**2)
         theta(a)=ACOS(top/bottom)
      END DO

      DO a=1,t
         i=da1(a)
         j=da2(a)
         k=da3(a)
         l=da4(a)
   
         colin=0
         ax=x(j)-x(i)
         ay=y(j)-y(i)
         az=z(j)-z(i)
         bx=x(k)-x(j)
         by=y(k)-y(j)
         bz=z(k)-z(j)
         cx=x(j)-x(l)
         cy=y(j)-y(l)
         cz=z(j)-z(l)
         dx=x(l)-x(k)
         dy=y(l)-y(k)
         dz=z(l)-z(k)
   
         IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) .AND. ay*bz.GT.(by*az-TINY) 
     1        .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
         IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) .AND. by*dz.GT.(dy*bz-TINY) 
     1        .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2

         IF (colin.EQ.1) THEN
            PRINT *,'Three sites',i,j,k,'are colinear'
              STOP
         ELSE IF (colin.EQ.2) THEN
            PRINT *,'Three sites',j,k,l,'are colinear'
              STOP
         ELSE
            lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
            mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
            qx(1)=x(i)+(lambda*bx)
            qy(1)=y(i)+(lambda*by)
            qz(1)=z(i)+(lambda*bz)
            qx(2)=x(j)
            qy(2)=y(j)
            qz(2)=z(j)
            qx(3)=x(l)+(mu*bx)
            qy(3)=y(l)+(mu*by)
            qz(3)=z(l)+(mu*bz)
            numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
            denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1          *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))
   
            IF ((numer/denom).LT. -1) numer=numer+1D-12
            IF ((numer/denom).GT.1) numer=numer-1D-12
  
            dphi(a)=ACOS(numer/denom)
         ENDIF
      END DO

      DO a=1,imp 
         e=ia1(a)
         f=ia2(a)
         g=ia3(a)
         h=ia4(a)    
         colin=0
         ax=x(f)-x(e)
         ay=y(f)-y(e)
         az=z(f)-z(e)
         bx=x(g)-x(f)
         by=y(g)-y(f)
         bz=z(g)-z(f)
         cx=x(f)-x(h)
         cy=y(f)-y(h)
         cz=z(f)-z(h)
         dx=x(h)-x(g)
         dy=y(h)-y(g)
         dz=z(h)-z(g)

         IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) .AND. ay*bz.GT.(by*az-TINY) 
     1      .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
         IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) .AND. by*dz.GT.(dy*bz-TINY) 
     1      .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2
         IF (colin.EQ.1) THEN
            PRINT *,'Three sites',e,f,g,' are colinear in improper ',a
            STOP
         ELSE IF (colin.EQ.2) THEN
            PRINT *,'Three sites',f,g,h,' are colinear in improper ',a
            STOP
         ELSE

            lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
            mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
            qx(1)=x(e)+(lambda*bx)
            qy(1)=y(e)+(lambda*by)
            qz(1)=z(e)+(lambda*bz)
            qx(2)=x(f)
            qy(2)=y(f)
            qz(2)=z(f)
            qx(3)=x(h)+(mu*bx)
            qy(3)=y(h)+(mu*by)
            qz(3)=z(h)+(mu*bz)
              numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
            denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1          *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))

            IF ((numer/denom).LT. -1) numer=numer+1D-12
            IF ((numer/denom).GT.1) numer=numer-1D-12

            iphi(a)=ACOS(numer/denom)
         ENDIF
      END DO

      RETURN

      END
