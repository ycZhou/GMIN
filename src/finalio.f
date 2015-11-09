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
      SUBROUTINE FINALIO
      USE commons
      use modamber
      use qmodule
      USE modcharmm
      IMPLICIT NONE

      INTEGER J1, J2, J3, NTYPEA
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB, RBCOORDS(NRBSITES*3), DCOORDS(3*NATOMS)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA

      DOUBLE PRECISION EPS2, RAD, HEIGHT
      COMMON /CAPS/ EPS2, RAD, HEIGHT

C dae rename RES to CHRES, TYPE to CHTYPE. Ugly fix to prevent conflict with amber

      INTEGER MAXSEG,MAXRES,MAXAIM
      PARAMETER (MAXSEG = 500)
      PARAMETER (MAXRES = 500)
      PARAMETER (MAXAIM = 510)
      CHARACTER(LEN=4) :: SEGID, CHRES, RESID, CHTYPE
      CHARACTER(LEN=6) :: CRMS
      COMMON /PSFC/ SEGID(MAXSEG), CHRES(MAXRES), RESID(MAXRES),
     &              CHTYPE(MAXAIM)
C     DOUBLE PRECISION CHX(MAXAIM),CHY(MAXAIM),CHZ(MAXAIM)
      CHARACTER(LEN=15) DBNAME(100)

C
      DBNAME(1)='dbase.1'
      DBNAME(2)='dbase.2'
      DBNAME(3)='dbase.3'
      DBNAME(4)='dbase.4'
      DBNAME(5)='dbase.5'
      DBNAME(6)='dbase.6'
      DBNAME(7)='dbase.7'
      DBNAME(8)='dbase.8'
      DBNAME(9)='dbase.9'
      DBNAME(10)='dbase.10'
      DBNAME(11)='dbase.11'
      DBNAME(12)='dbase.12'
      DBNAME(13)='dbase.13'
      DBNAME(14)='dbase.14'
      DBNAME(15)='dbase.15'
      DBNAME(16)='dbase.16'
      DBNAME(17)='dbase.17'
      DBNAME(18)='dbase.18'
      DBNAME(19)='dbase.19'
      DBNAME(20)='dbase.20'
      DBNAME(21)='dbase.21'
      DBNAME(22)='dbase.22'
      DBNAME(23)='dbase.23'
      DBNAME(24)='dbase.24'
      DBNAME(25)='dbase.25'
      DBNAME(26)='dbase.26'
      DBNAME(27)='dbase.27'
      DBNAME(28)='dbase.28'
      DBNAME(29)='dbase.29'
      DBNAME(30)='dbase.30'
      DBNAME(31)='dbase.31'
      DBNAME(32)='dbase.32'
      DBNAME(33)='dbase.33'
      DBNAME(34)='dbase.34'
      DBNAME(35)='dbase.35'
      DBNAME(36)='dbase.36'
      DBNAME(37)='dbase.37'
      DBNAME(38)='dbase.38'
      DBNAME(39)='dbase.39'
      DBNAME(40)='dbase.40'
      DBNAME(41)='dbase.41'
      DBNAME(42)='dbase.42'
      DBNAME(43)='dbase.43'
      DBNAME(44)='dbase.44'
      DBNAME(45)='dbase.45'
      DBNAME(46)='dbase.46'
      DBNAME(47)='dbase.47'
      DBNAME(48)='dbase.48'
      DBNAME(49)='dbase.49'
      DBNAME(50)='dbase.50'
      DBNAME(51)='dbase.51'
      DBNAME(52)='dbase.52'
      DBNAME(53)='dbase.53'
      DBNAME(54)='dbase.54'
      DBNAME(55)='dbase.55'
      DBNAME(56)='dbase.56'
      DBNAME(57)='dbase.57'
      DBNAME(58)='dbase.58'
      DBNAME(59)='dbase.59'
      DBNAME(60)='dbase.60'
      DBNAME(61)='dbase.61'
      DBNAME(62)='dbase.62'
      DBNAME(63)='dbase.63'
      DBNAME(64)='dbase.64'
      DBNAME(65)='dbase.65'
      DBNAME(66)='dbase.66'
      DBNAME(67)='dbase.67'
      DBNAME(68)='dbase.68'
      DBNAME(69)='dbase.69'
      DBNAME(70)='dbase.70'
      DBNAME(71)='dbase.71'
      DBNAME(72)='dbase.72'
      DBNAME(73)='dbase.73'
      DBNAME(74)='dbase.74'
      DBNAME(75)='dbase.75'
      DBNAME(76)='dbase.76'
      DBNAME(77)='dbase.77'
      DBNAME(78)='dbase.78'
      DBNAME(79)='dbase.79'
      DBNAME(80)='dbase.80'
      DBNAME(81)='dbase.81'
      DBNAME(82)='dbase.82'
      DBNAME(83)='dbase.83'
      DBNAME(84)='dbase.84'
      DBNAME(85)='dbase.85'
      DBNAME(86)='dbase.86'
      DBNAME(87)='dbase.87'
      DBNAME(88)='dbase.88'
      DBNAME(89)='dbase.89'
      DBNAME(90)='dbase.90'
      DBNAME(91)='dbase.91'
      DBNAME(92)='dbase.92'
      DBNAME(93)='dbase.93'
      DBNAME(94)='dbase.94'
      DBNAME(95)='dbase.95'
      DBNAME(96)='dbase.96'
      DBNAME(97)='dbase.97'
      DBNAME(98)='dbase.98'
      DBNAME(99)='dbase.99'
      DBNAME(100)='dbase.100'

      OPEN(UNIT=25,FILE='lowest',STATUS='UNKNOWN')
      DO J1=1,NSAVE
C        IF (FF(J1).EQ.0) STOP
         IF (RGCL2.OR.ARNO) THEN 
            WRITE(25,*) NATOMS+2
         ELSE
            WRITE(25,*) NATOMS
         ENDIF
         WRITE(25,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I2,'=',F20.10,' first found at step ',I8)
         IF (MSORIGT.OR.FRAUSIT) THEN
            WRITE(25,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
20          FORMAT('Si',3F20.10)
         ELSE IF (MSTRANST) THEN
            WRITE(25,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
         ELSE IF (RGCL2) THEN
            WRITE(25,'(A,F20.10)') 'Cl 0.0 0.0 ', 0.995D0
            WRITE(25,'(A,F20.10)') 'Cl 0.0 0.0 ',-0.995D0
            WRITE(25,60) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
60          FORMAT('AR ',3F20.10)
         ELSE IF (ARNO) THEN
            WRITE(25,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(25,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(25,65) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TOSI.OR.WELCH) THEN
            DO J2=1,NATOMS
               IF (ZSYM(J2).EQ.'PL') WRITE(25,'(A,3F20.10)') 'Na  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               IF (ZSYM(J2).EQ.'MI') WRITE(25,'(A,3F20.10)') 'Cl  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE IF (BLJCLUSTER) THEN
            DO J2=1,NATOMS
               IF (J2.LE.NTYPEA) THEN 
                  WRITE(25,'(A,3F20.10)') 'LA  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               ELSE
                  WRITE(25,'(A,3F20.10)') 'LB  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               ENDIF
            ENDDO
         ELSE IF (AMBER) THEN
            DO J2=1,NATOMS
               WRITE(25,'(A,3F20.10)') typech(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE IF (CHRMMT) THEN

C this writes 'lowest' in xyz (Xmakemol) format

            DO J2=1,NATOMS
               WRITE(25,'(A,1X,3F20.10)') ZSYM(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2=1,NATOMS
               DCOORDS(3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)
               DCOORDS(3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)
               DCOORDS(3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)
            ENDDO
            CALL CHARMMDUMP(DCOORDS,DBNAME(J1))

         ELSE IF (BLNT.AND.(.NOT.P46)) THEN
            DO J2=1,NATOMS
               WRITE(25,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            WRITE(25,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30          FORMAT('LA ',3F20.10)
         ENDIF
         IF ((NS.GT.0).AND.(.NOT.(WELCH.OR.TOSI))) THEN
            IF (MSORIGT.OR.FRAUSIT) THEN
               WRITE(25,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
40             FORMAT('Si',3F20.10)
            ELSE IF (MSTRANST) THEN
               WRITE(25,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
            ELSE
               WRITE(25,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
50             FORMAT('LB',3F20.10)
            ENDIF
         ENDIF
         IF (AMBER) CALL amberdump(J1,QMINP)
      ENDDO
      CLOSE(25)
 
      IF(CHRMMT.AND.RMST) THEN
        OPEN(UNIT=25,FILE='rmsbest',STATUS='UNKNOWN')
        DO J2=1,RMSSAVE
          WRITE(25,'(I6,F6.3,F15.5)')J2,RMSBEST(J2,1),RMSBEST(J2,2)
          WRITE(CRMS,'(I6)') J2
          DCOORDS(1:3*NATOMS)=RMSCOOR(J2,1:3*NATOMS)
          IF(RMSBEST(J2,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'rms.'//TRIM(ADJUSTL(CRMS)))
        ENDDO
        CLOSE(25)
      ENDIF

      IF (STOCKT) THEN
         OPEN(UNIT=26,FILE='stock.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS/2
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               WRITE(26,'(A4,3F18.10,A12,3F18.10)') 'LA ',
     &                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &                    'atom_vector',
     &                    SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &                    SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &                    COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))
            ENDDO
         ENDDO
         CLOSE(26)
      ELSE IF (TIP) THEN
         OPEN(UNIT=26,FILE='tip.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*3
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               CALL TIPIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1            QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                     RBCOORDS)
               WRITE(26,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ENDDO
         CLOSE(26)
      ELSE IF (CAPSID) THEN
         OPEN(UNIT=26,FILE='capsid.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*6
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               CALL CAPSIDIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                       RBCOORDS,RAD,HEIGHT)
               DO J3=1,5
                  WRITE(26,'(A4,3F20.10)') 'C1 ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
               WRITE(26,'(A4,3F20.10)') 'C4  ',RBCOORDS(16),RBCOORDS(17),RBCOORDS(18)
            ENDDO
         ENDDO
         CLOSE(26)
      ELSE IF (RIGID) THEN
         OPEN(UNIT=26,FILE='rigid.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               CALL RBIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),
     2                       QMINP(J1,3*(NATOMS/2+J2-1)+2),
     3                       QMINP(J1,3*(NATOMS/2+J2-1)+3),
     4                       RBCOORDS,NRBSITES,SITE)
               DO J3=1,NRBSITES
                  WRITE(26,'(A4,3F20.10)') 'LA ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(26)
      ENDIF

      RETURN
      END

      SUBROUTINE amberdump(J1,QMINP)
      USE commons
      USE modamber
      IMPLICIT NONE


      CHARACTER(LEN=25) coordfile
      CHARACTER(LEN=2) FNAME
      INTEGER J1
      DOUBLE PRECISION QMINP(NSAVE,3*NATOMS)

      IF (J1.LT.10) THEN
         WRITE (FNAME,'(I1)') J1
      ELSE
         WRITE (FNAME,'(I2)') J1
      ENDIF

      DO a=1,atoms
        x(a)=QMINP(J1,3*a-2)
        y(a)=QMINP(J1,3*a-1)
        z(a)=QMINP(J1,3*a)
      END DO

      coordfile='acoords.dump.'//FNAME

      OPEN (UNIT=4,IOSTAT=ios,FILE=coordfile,STATUS='UNKNOWN')

      DO a=1,atoms
        WRITE (UNIT=4,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F7.3,3X,F7.3,3X,F7.3)') label(a),typech(a),
     1        a,bondedto(a),x(a),y(a),z(a)
      ENDDO

      WRITE (UNIT=4,FMT='(A3)') 'end'
      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A4,7X,I2)') 'loop',rings

      DO a=1,rings
        WRITE (UNIT=4,FMT='(I3,4X,I3)') loopatom(2*a-1),loopatom(2*a)
      END DO

      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A7)') 'charges'

      DO a=1,atoms
        q(a)=q(a)/18.2223
        WRITE (UNIT=4,FMT='(I3,2X,F7.4)') a,q(a)
      END DO

      WRITE (UNIT=4,FMT='(A3)') 'end'

      RETURN

      END
C
C  SUBROUTINE to convert capsid CofM and DV coordinates to penatgons.
C
      SUBROUTINE CAPSIDIO(X1, Y1, Z1, L1, M1, N1,COORDS,RAD,HEIGHT)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), HEIGHT, C2A1,
     2                 M1, L1, N1, ALPHA1, RAD, CA1, S1, C3A1,
     3                 NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

      NUM1=-(1.0D0+SQRT(5.0D0))/4.0D0
      NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)/2.0D0
      NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)/2.0D0
      NUM4=(SQRT(5.0D0)-1.0D0)/4.0D0
      NUM5=-(1.0D0+SQRT(5.0D0))/4.0D0

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=RAD*CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=RAD*(-ALPHA1/2+ALPHA1**3/24)
         C3A1=RAD*(-0.5D0+ALPHA1**2/24.0D0)
         S1=RAD*(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=RAD*(CA1-1.0D0)/ALPHA1**2
         S1=RAD*SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) =     c2a1 - c3a1*l12 + x1
      COORDS(2) =     -(c3a1*l1*m1) - n1*s1 + y1
      COORDS(3) =     -(c3a1*l1*n1) + m1*s1 + z1
      COORDS(4) =     c2a1*num4 - c3a1*l1*(m1*num3 + l1*num4) + n1*num3*s1 + x1
      COORDS(5) =     c2a1*num3 - c3a1*m1*(m1*num3 + l1*num4) - n1*num4*s1 + y1
      COORDS(6) =     -(c3a1*n1*(m1*num3 + l1*num4)) - l1*num3*s1 + m1*num4*s1 + z1
      COORDS(7) =     c2a1*num1 - c3a1*l1*(l1*num1 + m1*num2) + n1*num2*s1 + x1
      COORDS(8) = c2a1*num2 - c3a1*m1*(l1*num1 + m1*num2) - n1*num5*s1 + y1
      COORDS(9) = -(c3a1*n1*(l1*num1 + m1*num2)) + m1*num1*s1 - l1*num2*s1 + z1
      COORDS(10) = c2a1*num1 + c3a1*l1*(-(l1*num1) + m1*num2) - n1*num2*s1 + x1
      COORDS(11) = -(c2a1*num2) + c3a1*m1*(-(l1*num1) + m1*num2) - n1*num5*s1 + y1
      COORDS(12) = -(c3a1*l1*n1*num1) + c3a1*m1*n1*num2 + m1*num1*s1 + l1*num2*s1 + z1
      COORDS(13) = c2a1*num4 + c3a1*l1*(m1*num3 - l1*num4) - n1*num3*s1 + x1
      COORDS(14) = -(c2a1*num3) + c3a1*m1*(m1*num3 - l1*num4) - n1*num4*s1 + y1
      COORDS(15) = c3a1*n1*(m1*num3 - l1*num4) + l1*num3*s1 + m1*num4*s1 + z1
C     COORDS(16)= (-(c3a1*l1*n1) - m1*s1 + 2*x1)/2.
C     COORDS(17)= -(c3a1*m1*n1)/2. + (l1*s1)/2. + y1
C     COORDS(18)= (c2a1 - c3a1*n12 + 2*z1)/2.
      COORDS(16)= -(c3a1*height*l1*n1) - height*m1*s1 + x1
      COORDS(17)= -(c3a1*height*m1*n1) + height*l1*s1 + y1
      COORDS(18)= c2a1*height - c3a1*height*n12 + z1

      RETURN
      END
C
C  Subroutine to convert rigid body CofM and DV coordinates to molecular sites.
C
      SUBROUTINE RBIO(X1, Y1, Z1, L1, M1, N1, COORDS, NRBSITES, SITE)
      IMPLICIT NONE
      INTEGER NRBSITES
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, SITE(NRBSITES,3),
     2                 M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12
      INTEGER J1

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF
   
      DO J1=1,NRBSITES
         COORDS(3*(J1-1)+1)=c2a1*SITE(J1,1) + s1*(n1*SITE(J1,2) - m1*SITE(J1,3)) - 
     1                      c3a1*l1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + X1
         COORDS(3*(J1-1)+2)=c2a1*SITE(J1,2) + s1*(-(n1*SITE(J1,1)) + l1*SITE(J1,3)) 
     1                    - c3a1*m1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Y1
         COORDS(3*(J1-1)+3)=s1*(m1*SITE(J1,1) - l1*SITE(J1,2)) + c2a1*SITE(J1,3) 
     1                    - c3a1*n1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Z1
      ENDDO 

      RETURN
      END
C
C  SUBROUTINE to convert TIP oxygen and DV coordinates to Cartesians.
C
      SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) = X1
      COORDS(2) = Y1
      COORDS(3) = Z1    
      COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
      COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
      COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
      COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

      RETURN
      END

