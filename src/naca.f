C  This file is part of libxfoil

C  libxfoil is free software: you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.

C  libxfoil is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.

C  Copyright (C) 2019 Daniel Prosser (this modified version of 
C  XFoil code)
C  Original copyright (C) 2000 Mark Drela (original XFoil code)

C===================================================================70

      SUBROUTINE NACA4(CAMB,XCAMB,THICK,XX,YT,YC,NSIDE,XB,YB,NB)
      REAL*8 CAMB, XCAMB, THICK
      REAL*8 XX(NSIDE), YT(NSIDE), YC(NSIDE)
      REAL*8 XB(2*NSIDE), YB(2*NSIDE)
      REAL*8 M
C
C---- TE point bunching parameter
      DATA AN / 1.5 /
C
      M = CAMB / 100.0
      P = XCAMB
      T = THICK / 100.0
C
      ANP = AN + 1.0
      DO 10 I=1, NSIDE
        FRAC = FLOAT(I-1)/FLOAT(NSIDE-1)
        IF(I.EQ.NSIDE) THEN
         XX(I) = 1.0
        ELSE
         XX(I) = 1.0 - ANP*FRAC*(1.0-FRAC)**AN - (1.0-FRAC)**ANP
        ENDIF
        YT(I) = ( 0.29690*SQRT(XX(I))
     &          - 0.12600*XX(I)
     &          - 0.35160*XX(I)**2
     &          + 0.28430*XX(I)**3
     &          - 0.10150*XX(I)**4) * T / 0.20
        IF(XX(I).LT.P) THEN
         YC(I) = M/P**2 * (2.0*P*XX(I) - XX(I)**2)
        ELSE
         YC(I) = M/(1.0-P)**2 * ((1.0-2.0*P) + 2.0*P*XX(I)-XX(I)**2)
        ENDIF
   10 CONTINUE
C
      IB = 0
      DO 20 I=NSIDE, 1, -1
        IB = IB + 1
        XB(IB) = XX(I)
        YB(IB) = YC(I) + YT(I)
   20 CONTINUE
      DO 30 I=2, NSIDE
        IB = IB + 1
        XB(IB) = XX(I)
        YB(IB) = YC(I) - YT(I)
   30 CONTINUE
      NB = IB
C
      RETURN
      END


      SUBROUTINE NACA5(IDES,XX,YT,YC,NSIDE,XB,YB,NB)
      REAL*8 XX(NSIDE), YT(NSIDE), YC(NSIDE)
      REAL*8 XB(2*NSIDE), YB(2*NSIDE)
      REAL*8 M
C
C---- TE point bunching parameter
      DATA AN / 1.5 /
C
      N5 =  IDES                                        / 10000
      N4 = (IDES - N5*10000                           ) / 1000
      N3 = (IDES - N5*10000 - N4*1000                 ) / 100
      N2 = (IDES - N5*10000 - N4*1000 - N3*100        ) / 10
      N1 = (IDES - N5*10000 - N4*1000 - N3*100 - N2*10)
C
      N543 = 100*N5 + 10*N4 + N3
C
      IF      (N543 .EQ. 210) THEN
cc     P = 0.05
       M = 0.0580
       C = 361.4
      ELSE IF (N543 .EQ. 220) THEN
cc     P = 0.10
       M = 0.1260
       C = 51.64
      ELSE IF (N543 .EQ. 230) THEN
cc     P = 0.15
       M = 0.2025
       C = 15.957
      ELSE IF (N543 .EQ. 240) THEN
cc     P = 0.20
       M = 0.2900
       C = 6.643
      ELSE IF (N543 .EQ. 250) THEN
cc     P = 0.25
       M = 0.3910
       C = 3.230
      ELSE
       WRITE(*,*) 'Illegal 5-digit designation'
       WRITE(*,*) 'First three digits must be 210, 220, ... 250'
       IDES = 0
       RETURN
      ENDIF
C
C
      T = FLOAT(N2*10 + N1) / 100.0
C
      ANP = AN + 1.0
      DO 10 I=1, NSIDE
        FRAC = FLOAT(I-1)/FLOAT(NSIDE-1)
        IF(I.EQ.NSIDE) THEN
         XX(I) = 1.0
        ELSE
         XX(I) = 1.0 - ANP*FRAC*(1.0-FRAC)**AN - (1.0-FRAC)**ANP
        ENDIF
C
        YT(I) = ( 0.29690*SQRT(XX(I))
     &          - 0.12600*XX(I)
     &          - 0.35160*XX(I)**2
     &          + 0.28430*XX(I)**3
     &          - 0.10150*XX(I)**4) * T / 0.20
        IF(XX(I).LT.M) THEN
         YC(I) = (C/6.0) * (XX(I)**3 - 3.0*M*XX(I)**2
     &                          + M*M*(3.0-M)*XX(I))
        ELSE
         YC(I) = (C/6.0) * M**3 * (1.0 - XX(I))
        ENDIF
   10 CONTINUE
C
      IB = 0
      DO 20 I=NSIDE, 1, -1
        IB = IB + 1
        XB(IB) = XX(I)
        YB(IB) = YC(I) + YT(I)
   20 CONTINUE
      DO 30 I=2, NSIDE
        IB = IB + 1
        XB(IB) = XX(I)
        YB(IB) = YC(I) - YT(I)
   30 CONTINUE
      NB = IB
C
      RETURN
      END
