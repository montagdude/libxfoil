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

C  Copyright (C) 2018 Daniel Prosser (this modified version of 
C  XFoil code)
C  Original copyright (C) 2000 Mark Drela (original XFoil code)

C===================================================================70
C
C     Sets actual Mach, Reynolds numbers
C     from unit-CL values and specified CLS
C     depending on MATYP,RETYP flags.
C
C===================================================================70
      SUBROUTINE MRCL(xfd,CLS,M_CLS,R_CLS)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd

      REAL*8 M_CLS
C
      CLA = MAX( CLS , 0.000001 )
C
C     DP mod: added SILENT_MODE option
      IF(xfd%RETYP.LT.1 .OR. xfd%RETYP.GT.3) THEN
        IF (.NOT. xfd%SILENT_MODE) THEN
          WRITE(*,*) 'MRCL:  Illegal Re(xfd%CL) dependence trigger.'
          WRITE(*,*) '       Setting fixed Re.'
        ENDIF
        xfd%RETYP = 1
      ENDIF
      IF(xfd%MATYP.LT.1 .OR. xfd%MATYP.GT.3) THEN
        IF (.NOT. xfd%SILENT_MODE) THEN
          WRITE(*,*) 'MRCL:  Illegal Mach(xfd%CL) dependence trigger.'
          WRITE(*,*) '       Setting fixed Mach.'
        ENDIF
        xfd%MATYP = 1
      ENDIF
C
C
      IF(xfd%MATYP.EQ.1) THEN
C
        xfd%MINF  = xfd%MINF1
        M_CLS = 0.
C
      ELSE IF(xfd%MATYP.EQ.2) THEN
C
        xfd%MINF  =  xfd%MINF1/SQRT(CLA)
        M_CLS = -0.5*xfd%MINF/CLA
C
      ELSE IF(xfd%MATYP.EQ.3) THEN
C
        xfd%MINF  = xfd%MINF1
        M_CLS = 0.
C
      ENDIF
C
C
      IF(xfd%RETYP.EQ.1) THEN
C
        xfd%REINF = xfd%REINF1
        R_CLS = 0.
C
      ELSE IF(xfd%RETYP.EQ.2) THEN
C
        xfd%REINF =  xfd%REINF1/SQRT(CLA)
        R_CLS = -0.5*xfd%REINF/CLA
C
      ELSE IF(xfd%RETYP.EQ.3) THEN
C
        xfd%REINF =  xfd%REINF1/CLA
        R_CLS = -xfd%REINF /CLA
C
      ENDIF
C
C
      IF(xfd%MINF .GE. 0.99) THEN
        IF (.NOT. xfd%SILENT_MODE) THEN
          WRITE(*,*)
          WRITE(*,*) 'MRCL: xfd%CL too low for chosen Mach(xfd%CL)
     &   dependence'
          WRITE(*,*) '      Aritificially limiting Mach to  0.99'
        ENDIF
        xfd%MINF = 0.99
        M_CLS = 0.
      ENDIF
C
      RRAT = 1.0
      IF(xfd%REINF1 .GT. 0.0) RRAT = xfd%REINF/xfd%REINF1
C
      IF(RRAT .GT. 100.0) THEN
        IF (.NOT. xfd%SILENT_MODE) THEN
          WRITE(*,*)
          WRITE(*,*) 'MRCL: xfd%CL too low for chosen Re(xfd%CL)
     &   dependence'
          WRITE(*,*) '      Aritificially limiting Re to ',xfd%REINF1
     &  *100.0
        ENDIF
        xfd%REINF = xfd%REINF1*100.0
        R_CLS = 0.
      ENDIF
C
      RETURN
      END ! MRCL

C===================================================================70
C
C     Sets the BL Newton system line number
C     corresponding to each BL station.
C
C===================================================================70
      SUBROUTINE IBLSYS(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      IV = 0
      DO 10 IS=1, 2
        DO 110 IBL=2, xfd%NBL(IS)
          IV = IV+1
          xfd%ISYS(IBL,IS) = IV
  110   CONTINUE
   10 CONTINUE
C
      xfd%NSYS = IV
      IF(xfd%NSYS.GT.2*IVX) STOP '*** IBLSYS: BL system array overflow. 
     &  ***'
C
      RETURN
      END

C===================================================================70
C 
C Set Karman-Tsien parameter TKLAM
C
C===================================================================70
      SUBROUTINE COMSET(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      BETA = SQRT(1.0 - xfd%MINF**2)
      BETA_MSQ = -0.5/BETA
C
      xfd%TKLAM   = xfd%MINF**2 / (1.0 + BETA)**2
      xfd%TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*xfd%TKLAM/ (1.0 + BETA) * BETA_MSQ
C
C---- set sonic Pressure coefficient and speed
      IF(xfd%MINF.EQ.0.0) THEN
       xfd%CPSTAR = -999.0
       xfd%QSTAR = 999.0
      ELSE
       xfd%CPSTAR = 2.0 / (xfd%GAMMA*xfd%MINF**2)
     &        * (( (1.0 + 0.5*xfd%GAMM1*xfd%MINF**2)
     &            /(1.0 + 0.5*xfd%GAMM1        ))**(xfd%GAMMA/xfd%GAMM1)
     &   - 1.0)
       xfd%QSTAR = xfd%QINF/xfd%MINF
     &       * SQRT( (1.0 + 0.5*xfd%GAMM1*xfd%MINF**2)
     &              /(1.0 + 0.5*xfd%GAMM1        ) )
      ENDIF
C
      RETURN
      END ! COMSET

C===================================================================70
C
C     Sets forced-transition BL coordinate locations.
C
C===================================================================70
      SUBROUTINE XIFSET(xfd,xbd,IS)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(xbl_data_type), intent(inout) :: xbd
C
      IF(xfd%XSTRIP(IS).GE.1.0) THEN
       xbd%XIFORC = xfd%XSSI(xfd%IBLTE(IS),IS)
       RETURN
      ENDIF
C
      CHX = xfd%XTE - xfd%XLE
      CHY = xfd%YTE - xfd%YLE
      CHSQ = CHX**2 + CHY**2
C
C---- calculate chord-based x/c, y/c
      DO 10 I=1, xfd%N
        xfd%W1(I) = ((xfd%X(I)-xfd%XLE)*CHX + (xfd%Y(I)-xfd%YLE)*CHY) /
     &   CHSQ
        xfd%W2(I) = ((xfd%Y(I)-xfd%YLE)*CHX - (xfd%X(I)-xfd%XLE)*CHY) /
     &   CHSQ
 10   CONTINUE
C
      CALL SPLIND(xfd%W1,xfd%W3,xfd%S,xfd%N,-999.0,-999.0)
      CALL SPLIND(xfd%W2,xfd%W4,xfd%S,xfd%N,-999.0,-999.0)
C
      IF(IS.EQ.1) THEN
C
C----- set approximate arc length of forced transition point for SINVRT
       STR = xfd%SLE + (xfd%S(1)-xfd%SLE)*xfd%XSTRIP(IS)
C
C----- calculate actual arc length
       CALL SINVRT(STR,xfd%XSTRIP(IS),xfd%W1,xfd%W3,xfd%S,xfd%N
     &  ,xfd%SILENT_MODE)
C
C----- set BL coordinate value
       xbd%XIFORC = MIN( (xfd%SST - STR) , xfd%XSSI(xfd%IBLTE(IS),IS) )
C
      ELSE
C----- same for bottom side
C
       STR = xfd%SLE + (xfd%S(xfd%N)-xfd%SLE)*xfd%XSTRIP(IS)
       CALL SINVRT(STR,xfd%XSTRIP(IS),xfd%W1,xfd%W3,xfd%S,xfd%N
     &  ,xfd%SILENT_MODE)
       xbd%XIFORC = MIN( (STR - xfd%SST) , xfd%XSSI(xfd%IBLTE(IS),IS) )
C
      ENDIF
C
      IF(xbd%XIFORC .LT. 0.0) THEN
C      DP mod: added SILENT_MODE option
       IF (.NOT. xfd%SILENT_MODE) WRITE(*,1000) IS
 1000  FORMAT(/' ***  Stagnation point is past trip on side',I2,'  ***')
       xbd%XIFORC = xfd%XSSI(xfd%IBLTE(IS),IS)
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Set BL primary "2" variables from parameter list
C
C===================================================================70
      SUBROUTINE BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(xbl_data_type), intent(inout) :: xbd
C
      xbd%X2 = XSI
      xbd%AMPL2 = AMI
      xbd%S2  = CTI
      xbd%T2  = THI
      xbd%D2  = DSI - DSWAKI
      xbd%DW2 = DSWAKI
C
      xbd%U2 = UEI*(1.0-xbd%TKBL) / (1.0 - xbd%TKBL*(UEI/xbd%QINFBL)**2)
     &  
      xbd%U2_UEI = (1.0 + xbd%TKBL*(2.0*xbd%U2*UEI/xbd%QINFBL**2 - 1.0))
     &  
     &       / (1.0 - xbd%TKBL*(UEI/xbd%QINFBL)**2)
      xbd%U2_MS  = (xbd%U2*(UEI/xbd%QINFBL)**2  -  UEI)*xbd%TKBL_MS
     &                    / (1.0 - xbd%TKBL*(UEI/xbd%QINFBL)**2)
C
      RETURN
      END ! BLPRV

C===================================================================70
C===================================================================70
      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )

      REAL*8 MSQ
C
C---- calculate kinematic shape parameter (assuming air)
C     (from Whitfield )
      HK     =    (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =     1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK)/(1.0 + 0.113*MSQ)
C
      RETURN
      END
C===================================================================70
C
C     Calculates turbulence-independent secondary "2" 
C     variables from the primary "2" variables.
C
C===================================================================70
      SUBROUTINE BLKIN(xbd)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(xbl_data_type), intent(inout) :: xbd

C     DP mod: set explicitly (otherwise initialized as 0 here)
      HVRAT = 0.0
C
C---- set edge Mach number ** 2
      xbd%M2    = xbd%U2*xbd%U2*xbd%HSTINV / (xbd%GM1BL*(1.0 - 0.5
     &  *xbd%U2*xbd%U2*xbd%HSTINV))
      TR2   = 1.0 + 0.5*xbd%GM1BL*xbd%M2
      xbd%M2_U2 = 2.0*xbd%M2*TR2/xbd%U2
      xbd%M2_MS = xbd%U2*xbd%U2*TR2    / (xbd%GM1BL*(1.0 - 0.5*xbd%U2
     &  *xbd%U2*xbd%HSTINV))
     &      * xbd%HSTINV_MS
C
C---- set edge static density (isentropic relation)
      xbd%R2    = xbd%RSTBL   *TR2**(-1.0/xbd%GM1BL)
      xbd%R2_U2 = -xbd%R2/TR2 * 0.5*xbd%M2_U2
      xbd%R2_MS = -xbd%R2/TR2 * 0.5*xbd%M2_MS
     &      + xbd%RSTBL_MS*TR2**(-1.0/xbd%GM1BL)
C
C---- set shape parameter
      xbd%H2    =  xbd%D2/xbd%T2
      xbd%H2_D2 = 1.0/xbd%T2
      xbd%H2_T2 = -xbd%H2/xbd%T2
C
C---- set edge static/stagnation enthalpy
      HERAT = 1.0 - 0.5*xbd%U2*xbd%U2*xbd%HSTINV
      HE_U2 =     -        xbd%U2*xbd%HSTINV
      HE_MS =     - 0.5*xbd%U2*xbd%U2*xbd%HSTINV_MS
C
C---- set molecular viscosity
      xbd%V2 = SQRT((HERAT)**3) * (1.0+HVRAT)/(HERAT+HVRAT)/xbd%REYBL
      V2_HE = xbd%V2*(1.5/HERAT - 1.0/(HERAT+HVRAT))
C
      xbd%V2_U2 =                        V2_HE*HE_U2
      xbd%V2_MS = -xbd%V2/xbd%REYBL * xbd%REYBL_MS + V2_HE*HE_MS
      xbd%V2_RE = -xbd%V2/xbd%REYBL * xbd%REYBL_RE
C
C---- set kinematic shape parameter
      CALL HKIN( xbd%H2, xbd%M2, xbd%HK2, HK2_H2, HK2_M2 )
C
      xbd%HK2_U2 =                HK2_M2*xbd%M2_U2
      xbd%HK2_T2 = HK2_H2*xbd%H2_T2
      xbd%HK2_D2 = HK2_H2*xbd%H2_D2
      xbd%HK2_MS =                HK2_M2*xbd%M2_MS
C
C---- set momentum thickness Reynolds number
      xbd%RT2    = xbd%R2*xbd%U2*xbd%T2/xbd%V2
      xbd%RT2_U2 = xbd%RT2*(1.0/xbd%U2 + xbd%R2_U2/xbd%R2 - xbd%V2_U2
     &  /xbd%V2)
      xbd%RT2_T2 = xbd%RT2/xbd%T2
      xbd%RT2_MS = xbd%RT2*(         xbd%R2_MS/xbd%R2 - xbd%V2_MS/xbd%V2
     &  )
      xbd%RT2_RE = xbd%RT2*(                  - xbd%V2_RE/xbd%V2)
C
      RETURN
      END ! BLKIN

C===================================================================70
C
C     Amplification rate routine for envelope e^n method.
C     Reference:
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and
C                Low Reynolds Number Airfoils",
C                AIAA Journal, Oct. 1987.
C
C     NEW VERSION.   March 1991       (latest bug fix  July 93)
C          - m(H) correlation made more accurate up to H=20
C          - for H > 5, non-similar profiles are used 
C            instead of Falkner-Skan profiles.  These 
C            non-similar profiles have smaller reverse 
C            velocities, are more representative of typical 
C            separation bubble profiles.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is
C            calculated by integrating AX (= dN/dx) with
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C
C===================================================================70
      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
      IMPLICIT REAL (A-H,M,O-Z)
ccc   DATA DGR / 0.04 /
      DATA DGR / 0.08 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) - H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRCRIT = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRCRIT-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRCRIT-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C----- Amplification envelope slope correlation for Falkner-Skan
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C----- new m(H) correlation    1 March 91
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2
       AF_HK = AF_HMI*HMI_HK
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      RETURN
      END ! DAMPL

C===================================================================70 
C
C     Amplification rate routine for modified envelope e^n method.
C     Reference: 
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and 
C                Low Reynolds Number Airfoils", 
C                AIAA Journal, Oct. 1987.
C
C     NEWER VERSION.   Nov 1996
C          - Amplification rate changes to the Orr-Sommerfeld
C              maximum ai(H,Rt) function for H > 4 .
C          - This implicitly assumes that the frequency range
C              (around w = 0.09 Ue/theta) which experiences this 
C              maximum amplification rate contains the currently
C              most-amplified frequency.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is 
C            calculated by integrating AX (= dN/dx) with 
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C
C===================================================================70 
      SUBROUTINE DAMPL2( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
      IMPLICIT REAL (A-H,M,O-Z)
      DATA DGR / 0.08 /
      DATA HK1, HK2 / 3.5, 4.0 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) -- H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRC = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRC-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRC-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C
C----- set envelope amplification rate with respect to Rtheta
C-       DADR = d(N)/d(Rtheta) = f(H)
C
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C
C----- set conversion factor from d/d(Rtheta) to d/dx
C-       AF = Theta d(Rtheta)/dx = f(H)
C
       BRG = -20.0*HMI
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3 + 0.1*EXP(BRG)
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2 - 2.0*EXP(BRG)
       AF_HK = AF_HMI*HMI_HK
C
C
C----- set amplification rate with respect to x, 
C-     with RFAC shutting off amplification when below Rcrit
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      IF(HK .LT. HK1) RETURN
C
C---- non-envelope max-amplification correction for separated profiles
C
      HNORM = (HK - HK1) / (HK2 - HK1)
      HN_HK =       1.0  / (HK2 - HK1)
C
C---- set blending fraction HFAC = 0..1 over HK1 < HK < HK2
      IF(HNORM .GE. 1.0) THEN
       HFAC = 1.0
       HF_HK = 0.
      ELSE
       HFAC  =  3.0*HNORM**2 - 2.0*HNORM**3
       HF_HK = (6.0*HNORM    - 6.0*HNORM**2)*HN_HK
      ENDIF
C
C---- "normal" envelope amplification rate AX1
      AX1    = AX
      AX1_HK = AX_HK
      AX1_TH = AX_TH
      AX1_RT = AX_RT
C
C---- set modified amplification rate AX2
      GR0 = 0.30 + 0.35 * EXP(-0.15*(HK-5.0))
      GR0_HK =   - 0.35 * EXP(-0.15*(HK-5.0)) * 0.15
C
      TNR = TANH(1.2*(GR - GR0))
      TNR_RT =  (1.0 - TNR**2)*1.2*GR_RT
      TNR_HK = -(1.0 - TNR**2)*1.2*GR0_HK
C
      AX2    = (0.086*TNR    -     0.25/(HK-1.0)**1.5) / TH
      AX2_HK = (0.086*TNR_HK + 1.5*0.25/(HK-1.0)**2.5) / TH
      AX2_RT = (0.086*TNR_RT                         ) / TH
      AX2_TH = -AX2/TH
C
      IF(AX2 .LT. 0.0) THEN
       AX2    = 0.0
       AX2_HK = 0.
       AX2_RT = 0.
       AX2_TH = 0.
      ENDIF
C
C---- blend the two amplification rates
      AX    = HFAC*AX2    + (1.0 - HFAC)*AX1
      AX_HK = HFAC*AX2_HK + (1.0 - HFAC)*AX1_HK + HF_HK*(AX2-AX1)
      AX_RT = HFAC*AX2_RT + (1.0 - HFAC)*AX1_RT
      AX_TH = HFAC*AX2_TH + (1.0 - HFAC)*AX1_TH
C
      RETURN
      END ! DAMPL2

C===================================================================70
C
C     Returns average amplification AX over interval 1..2
C
C===================================================================70
      SUBROUTINE AXSET( HK1,    T1,    RT1,    A1,
     &                  HK2,    T2,    RT2,    A2,  ACRIT, IDAMPV,
     &           AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &               AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
cC==========================
cC---- 1st-order -- based on "1" quantities only
c      CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
c      AX2_HK2 = 0.0
c      AX2_T2  = 0.0
c      AX2_RT2 = 0.0
cC
c      AX1_A1 = 0.0
c      AX2_A2 = 0.0
cC
c      AX     = AX1
c      AX_AX1 = 1.0
c      AX_AX2 = 0.0
cC
c      ARG = MIN( 20.0*(ACRIT-A1) , 20.0 )
c      EXN    = EXP(-ARG)
c      EXN_A1 = 20.0*EXN
c      EXN_A2 = 0.
cC
c      DAX    = EXN   * 0.0004/T1
c      DAX_A1 = EXN_A1* 0.0004/T1
c      DAX_A2 = 0.
c      DAX_T1 = -DAX/T1
c      DAX_T2 = 0.
C
C==========================
C---- 2nd-order
      IF(IDAMPV.EQ.0) THEN
       CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ELSE
       CALL DAMPL2( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL2( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ENDIF
C
CC---- simple-average version
C      AXA = 0.5*(AX1 + AX2)
C      IF(AXA .LE. 0.0) THEN
C       AXA = 0.0
C       AXA_AX1 = 0.0
C       AXA_AX2 = 0.0
C      ELSE
C       AXA_AX1 = 0.5
C       AXA_AX2 = 0.5
C      ENDIF
C
C---- rms-average version (seems a little better on coarse grids)
      AXSQ = 0.5*(AX1**2 + AX2**2)
      IF(AXSQ .LE. 0.0) THEN
       AXA = 0.0
       AXA_AX1 = 0.0
       AXA_AX2 = 0.0
      ELSE
       AXA = SQRT(AXSQ)
       AXA_AX1 = 0.5*AX1/AXA
       AXA_AX2 = 0.5*AX2/AXA
      ENDIF
C
C----- small additional term to ensure  dN/dx > 0  near  N = Ncrit
       ARG = MIN( 20.0*(ACRIT-0.5*(A1+A2)) , 20.0 )
       IF(ARG.LE.0.0) THEN
        EXN    = 1.0
CC      EXN_AC = 0.
        EXN_A1 = 0.
        EXN_A2 = 0.
       ELSE
        EXN    = EXP(-ARG)
CC      EXN_AC = -20.0    *EXN
        EXN_A1 =  20.0*0.5*EXN
        EXN_A2 =  20.0*0.5*EXN
       ENDIF
C
       DAX    = EXN    * 0.002/(T1+T2)
CC     DAX_AC = EXN_AC * 0.002/(T1+T2)
       DAX_A1 = EXN_A1 * 0.002/(T1+T2)
       DAX_A2 = EXN_A2 * 0.002/(T1+T2)
       DAX_T1 = -DAX/(T1+T2)
       DAX_T2 = -DAX/(T1+T2)
C
c
c        DAX    = 0.
c        DAX_A1 = 0.
c        DAX_A2 = 0.
c        DAX_AC = 0.
c        DAX_T1 = 0.
c        DAX_T2 = 0.
C==========================
C
      AX     = AXA             + DAX
C
      AX_HK1 = AXA_AX1*AX1_HK1
      AX_T1  = AXA_AX1*AX1_T1  + DAX_T1
      AX_RT1 = AXA_AX1*AX1_RT1
      AX_A1  =                   DAX_A1
C
      AX_HK2 = AXA_AX2*AX2_HK2
      AX_T2  = AXA_AX2*AX2_T2  + DAX_T2
      AX_RT2 = AXA_AX2*AX2_RT2
      AX_A2  =                   DAX_A2
C
      RETURN
      END

C===================================================================70
      SUBROUTINE TRCHEK2(xfd,xbd)
C
C     New second-order version:  December 1994.
C
C     Checks if transition occurs in the current interval X1..X2.
C     If transition occurs, then set transition location XT, and 
C     its sensitivities to "1" and "2" variables.  If no transition, 
C     set amplification AMPL2.
C
C
C     Solves the implicit amplification equation for N2:
C
C       N2 - N1     N'(XT,NT) + N'(X1,N1)
C       -------  =  ---------------------
C       X2 - X1               2
C
C     In effect, a 2-point central difference is used between
C     X1..X2 (no transition), or X1..XT (transition).  The switch
C     is done by defining XT,NT in the equation above depending
C     on whether N2 exceeds Ncrit.
C
C  If N2<Ncrit:  NT=N2    , XT=X2                  (no transition)
C
C  If N2>Ncrit:  NT=Ncrit , XT=(Ncrit-N1)/(N2-N1)  (transition)
C
C
C
      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(xbl_data_type), intent(inout) :: xbd
      DATA DAEPS / 5.0E-5 /
CCC   DATA DAEPS / 1.0D-12 /
C
C---- save variables and sensitivities at IBL ("2") for future restoration
C     DP mod: to remove need for EQUIVALENCE and COM2
      call store_c2sav(xbd)
C      DO 5 ICOM=1, NCOM
C        C2SAV(ICOM) = COM2(ICOM)
C    5 CONTINUE
C
C---- calculate average amplification rate AX over X1..X2 interval
      CALL AXSET( xbd%HK1,    xbd%T1,    xbd%RT1, xbd%AMPL1,
     &            xbd%HK2,    xbd%T2,    xbd%RT2, xbd%AMPL2,  xbd%AMCRIT
     &  , xbd%IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
C---- set initial guess for iterate N2 (AMPL2) at X2
      xbd%AMPL2 = xbd%AMPL1 + AX*(xbd%X2-xbd%X1)
C
C---- solve implicit system for amplification AMPL2
      DO 100 ITAM=1, 30
C
C---- define weighting factors WF1,WF2 for defining "T" quantities from 1,2
C
      IF(xbd%AMPL2 .LE. xbd%AMCRIT) THEN
C------ there is no transition yet,  "T" is the same as "2"
        AMPLT    = xbd%AMPL2
        AMPLT_A2 = 1.0
        SFA    = 1.0
        SFA_A1 = 0.
        SFA_A2 = 0.
      ELSE
C------ there is transition in X1..X2, "T" is set from N1, N2
        AMPLT    = xbd%AMCRIT
        AMPLT_A2 = 0.
        SFA    = (AMPLT - xbd%AMPL1)/(xbd%AMPL2-xbd%AMPL1)
        SFA_A1 = ( SFA  - 1.0  )/(xbd%AMPL2-xbd%AMPL1)
        SFA_A2 = (      - SFA  )/(xbd%AMPL2-xbd%AMPL1)
      ENDIF
C
      IF(xbd%XIFORC.LT.xbd%X2) THEN
        SFX    = (xbd%XIFORC - xbd%X1 )/(xbd%X2-xbd%X1)
        SFX_X1 = (SFX    - 1.0)/(xbd%X2-xbd%X1)
        SFX_X2 = (       - SFX)/(xbd%X2-xbd%X1)
        SFX_XF =  1.0          /(xbd%X2-xbd%X1)
      ELSE
        SFX    = 1.0
        SFX_X1 = 0.
        SFX_X2 = 0.
        SFX_XF = 0.
      ENDIF
C
C---- set weighting factor from free or forced transition
      IF(SFA.LT.SFX) THEN
        WF2    = SFA
        WF2_A1 = SFA_A1
        WF2_A2 = SFA_A2
        WF2_X1 = 0.
        WF2_X2 = 0.
        WF2_XF = 0.
      ELSE
        WF2    = SFX
        WF2_A1 = 0.
        WF2_A2 = 0.
        WF2_X1 = SFX_X1
        WF2_X2 = SFX_X2
        WF2_XF = SFX_XF
      ENDIF
C
C
C=====================
CC---- 1st-order (based on "1" quantites only, for testing)
C      WF2    = 0.0
C      WF2_A1 = 0.0
C      WF2_A2 = 0.0
C      WF2_X1 = 0.0
C      WF2_X2 = 0.0
C      WF2_XF = 0.0
C=====================
C
      WF1    = 1.0 - WF2
      WF1_A1 =     - WF2_A1
      WF1_A2 =     - WF2_A2
      WF1_X1 =     - WF2_X1
      WF1_X2 =     - WF2_X2
      WF1_XF =     - WF2_XF
C
C---- interpolate BL variables to XT
      xbd%XT    = xbd%X1*WF1    + xbd%X2*WF2
      TT    = xbd%T1*WF1    + xbd%T2*WF2
      DT    = xbd%D1*WF1    + xbd%D2*WF2
      UT    = xbd%U1*WF1    + xbd%U2*WF2
C
      XT_A2 = xbd%X1*WF1_A2 + xbd%X2*WF2_A2
      TT_A2 = xbd%T1*WF1_A2 + xbd%T2*WF2_A2
      DT_A2 = xbd%D1*WF1_A2 + xbd%D2*WF2_A2
      UT_A2 = xbd%U1*WF1_A2 + xbd%U2*WF2_A2
C
C---- temporarily set "2" variables from "T" for BLKIN
      xbd%X2 = xbd%XT
      xbd%T2 = TT
      xbd%D2 = DT
      xbd%U2 = UT
C
C---- calculate laminar secondary "T" variables HKT, RTT
      CALL BLKIN(xbd)
C
      HKT    = xbd%HK2
      HKT_TT = xbd%HK2_T2
      HKT_DT = xbd%HK2_D2
      HKT_UT = xbd%HK2_U2
      HKT_MS = xbd%HK2_MS
C
      RTT    = xbd%RT2
      RTT_TT = xbd%RT2_T2
      RTT_UT = xbd%RT2_U2
      RTT_MS = xbd%RT2_MS
      RTT_RE = xbd%RT2_RE
C
C---- restore clobbered "2" variables, except for AMPL2
      AMSAVE = xbd%AMPL2
C     DP mod: to remove need for EQUIVALENCE and COM2
      call from_c2sav(xbd)
C      DO 8 ICOM=1, NCOM
C        COM2(ICOM) = C2SAV(ICOM)
C 8    CONTINUE
      xbd%AMPL2 = AMSAVE
C
C---- calculate amplification rate AX over current X1-XT interval
      CALL AXSET( xbd%HK1,    xbd%T1,    xbd%RT1, xbd%AMPL1,
     &            HKT,    TT,    RTT, AMPLT,  xbd%AMCRIT, xbd%IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HKT, AX_TT, AX_RTT, AX_AT )
C
C---- punch out early if there is no amplification here
      IF(AX .LE. 0.0) GO TO 101
C
C---- set sensitivity of AX(A2)
      AX_A2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
     &      +  AX_AT                                 *AMPLT_A2
C
C---- residual for implicit AMPL2 definition (amplification equation)
      RES    = xbd%AMPL2 - xbd%AMPL1 - AX   *(xbd%X2-xbd%X1) 
      RES_A2 = 1.0           - AX_A2*(xbd%X2-xbd%X1)
C
      DA2 = -RES/RES_A2
C
      RLX = 1.0
      DXT = XT_A2*DA2
C
      IF(RLX*ABS(DXT/(xbd%X2-xbd%X1)) .GT. 0.05) RLX = 0.05*ABS((xbd%X2
     &  -xbd%X1)/DXT)
      IF(RLX*ABS(DA2)         .GT. 1.0 ) RLX = 1.0 *ABS(   1.0 /DA2)
C
C---- check if converged
      IF(ABS(DA2) .LT. DAEPS) GO TO 101
C
      IF((xbd%AMPL2.GT.xbd%AMCRIT .AND. xbd%AMPL2+RLX*DA2.LT.xbd%AMCRIT)
     &  .OR.
     &   (xbd%AMPL2.LT.xbd%AMCRIT .AND. xbd%AMPL2+RLX*DA2.GT.xbd%AMCRIT)
     &      ) THEN
C------ limited Newton step so AMPL2 doesn't step across AMCRIT either way
        xbd%AMPL2 = xbd%AMCRIT
      ELSE
C------ regular Newton step
        xbd%AMPL2 = xbd%AMPL2 + RLX*DA2
      ENDIF
C
 100  CONTINUE
C     DP mod: check for infinite loop condition
      IF (ISNAN(xbd%AMPL1) .OR. ISNAN(AMPLT) .OR. ISNAN(xbd%AMPL2)) THEN
     &  
        xfd%XFOIL_FAIL = .TRUE.
        RETURN
      ENDIF

C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE) THEN
        WRITE(*,*) 'TRCHEK2: N2 convergence failed.'
        WRITE(*,6700) xbd%X1, xbd%XT, xbd%X2, xbd%AMPL1, AMPLT,
     &   xbd%AMPL2, AX, DA2
      ENDIF
 6700 FORMAT(1X,'x:', 3F9.5,'  N:',3F7.3,'  Nx:',F8.3,'   dN:',E10.3)
C
 101  CONTINUE
C
C
C---- test for free or forced transition
      xbd%TRFREE = xbd%AMPL2 .GE. xbd%AMCRIT
      xbd%TRFORC = xbd%XIFORC.GT.xbd%X1 .AND. xbd%XIFORC.LE.xbd%X2
C
C---- set transition interval flag
      xbd%TRAN = xbd%TRFORC .OR. xbd%TRFREE
C
      IF(.NOT.xbd%TRAN) RETURN
C
C---- resolve if both forced and free transition
      IF(xbd%TRFREE .AND. xbd%TRFORC) THEN
       xbd%TRFORC = xbd%XIFORC .LT. xbd%XT
       xbd%TRFREE = xbd%XIFORC .GE. xbd%XT
      ENDIF
C
      IF(xbd%TRFORC) THEN
C----- if forced transition, then XT is prescribed,
C-     no sense calculating the sensitivities, since we know them...
       xbd%XT = xbd%XIFORC
       xbd%XT_A1 = 0.
       xbd%XT_X1 = 0.
       xbd%XT_T1 = 0.
       xbd%XT_D1 = 0.
       xbd%XT_U1 = 0.
       xbd%XT_X2 = 0.
       xbd%XT_T2 = 0.
       xbd%XT_D2 = 0.
       xbd%XT_U2 = 0.
       xbd%XT_MS = 0.
       xbd%XT_RE = 0.
       xbd%XT_XF = 1.0
       RETURN
      ENDIF
C
C---- free transition ... set sensitivities of XT
C
C---- XT( X1 X2 A1 A2 XF ),  TT( T1 T2 A1 A2 X1 X2 XF),   DT( ...
CC    XT    = X1*WF1    + X2*WF2
CC    TT    = T1*WF1    + T2*WF2
CC    DT    = D1*WF1    + D2*WF2
CC    UT    = U1*WF1    + U2*WF2
C
      xbd%XT_X1 =    WF1
      TT_T1 =    WF1
      DT_D1 =    WF1
      UT_U1 =    WF1
C
      xbd%XT_X2 =                WF2
      TT_T2 =                WF2
      DT_D2 =                WF2
      UT_U2 =                WF2
C
      xbd%XT_A1 = xbd%X1*WF1_A1 + xbd%X2*WF2_A1
      TT_A1 = xbd%T1*WF1_A1 + xbd%T2*WF2_A1
      DT_A1 = xbd%D1*WF1_A1 + xbd%D2*WF2_A1
      UT_A1 = xbd%U1*WF1_A1 + xbd%U2*WF2_A1
C
CC    XT_A2 = X1*WF1_A2 + X2*WF2_A2
CC    TT_A2 = T1*WF1_A2 + T2*WF2_A2
CC    DT_A2 = D1*WF1_A2 + D2*WF2_A2
CC    UT_A2 = U1*WF1_A2 + U2*WF2_A2
C
      xbd%XT_X1 = xbd%X1*WF1_X1 + xbd%X2*WF2_X1 + xbd%XT_X1
      TT_X1 = xbd%T1*WF1_X1 + xbd%T2*WF2_X1
      DT_X1 = xbd%D1*WF1_X1 + xbd%D2*WF2_X1
      UT_X1 = xbd%U1*WF1_X1 + xbd%U2*WF2_X1
C
      xbd%XT_X2 = xbd%X1*WF1_X2 + xbd%X2*WF2_X2 + xbd%XT_X2
      TT_X2 = xbd%T1*WF1_X2 + xbd%T2*WF2_X2
      DT_X2 = xbd%D1*WF1_X2 + xbd%D2*WF2_X2
      UT_X2 = xbd%U1*WF1_X2 + xbd%U2*WF2_X2
C
      xbd%XT_XF = xbd%X1*WF1_XF + xbd%X2*WF2_XF
      TT_XF = xbd%T1*WF1_XF + xbd%T2*WF2_XF
      DT_XF = xbd%D1*WF1_XF + xbd%D2*WF2_XF
      UT_XF = xbd%U1*WF1_XF + xbd%U2*WF2_XF
C
C---- at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT )
C
C---- set sensitivities of AX( T1 D1 U1 A1 T2 D2 U2 A2 MS RE )
      AX_T1 =  AX_HK1*xbd%HK1_T1 + AX_T1 + AX_RT1*xbd%RT1_T1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T1
      AX_D1 =  AX_HK1*xbd%HK1_D1
     &      + (AX_HKT*HKT_DT                        )*DT_D1
      AX_U1 =  AX_HK1*xbd%HK1_U1         + AX_RT1*xbd%RT1_U1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U1
      AX_A1 =  AX_A1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A1
     &      + (AX_HKT*HKT_DT                        )*DT_A1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A1
      AX_X1 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X1
     &      + (AX_HKT*HKT_DT                        )*DT_X1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X1
C
      AX_T2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T2
      AX_D2 = (AX_HKT*HKT_DT                        )*DT_D2
      AX_U2 = (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U2
      AX_A2 =  AX_AT                                 *AMPLT_A2
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
      AX_X2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X2
     &      + (AX_HKT*HKT_DT                        )*DT_X2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X2
C
      AX_XF = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_XF
     &      + (AX_HKT*HKT_DT                        )*DT_XF
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_XF
C
      AX_MS =  AX_HKT*HKT_MS         + AX_RTT*RTT_MS
     &      +  AX_HK1*xbd%HK1_MS         + AX_RT1*xbd%RT1_MS
      AX_RE =                          AX_RTT*RTT_RE
     &                               + AX_RT1*xbd%RT1_RE
C
C
C---- set sensitivities of residual RES
CCC   RES  = AMPL2 - AMPL1 - AX*(X2-X1)
      Z_AX =               -    (xbd%X2-xbd%X1)
C
      Z_A1 = Z_AX*AX_A1 - 1.0
      Z_T1 = Z_AX*AX_T1
      Z_D1 = Z_AX*AX_D1
      Z_U1 = Z_AX*AX_U1
      Z_X1 = Z_AX*AX_X1 + AX
C
      Z_A2 = Z_AX*AX_A2 + 1.0
      Z_T2 = Z_AX*AX_T2
      Z_D2 = Z_AX*AX_D2
      Z_U2 = Z_AX*AX_U2
      Z_X2 = Z_AX*AX_X2 - AX
C
      Z_XF = Z_AX*AX_XF
      Z_MS = Z_AX*AX_MS
      Z_RE = Z_AX*AX_RE
C
C---- set sensitivities of XT, with RES being stationary for A2 constraint
      xbd%XT_A1 = xbd%XT_A1 - (XT_A2/Z_A2)*Z_A1
      xbd%XT_T1 =       - (XT_A2/Z_A2)*Z_T1
      xbd%XT_D1 =       - (XT_A2/Z_A2)*Z_D1
      xbd%XT_U1 =       - (XT_A2/Z_A2)*Z_U1
      xbd%XT_X1 = xbd%XT_X1 - (XT_A2/Z_A2)*Z_X1
      xbd%XT_T2 =       - (XT_A2/Z_A2)*Z_T2
      xbd%XT_D2 =       - (XT_A2/Z_A2)*Z_D2
      xbd%XT_U2 =       - (XT_A2/Z_A2)*Z_U2
      xbd%XT_X2 = xbd%XT_X2 - (XT_A2/Z_A2)*Z_X2
      xbd%XT_MS =       - (XT_A2/Z_A2)*Z_MS
      xbd%XT_RE =       - (XT_A2/Z_A2)*Z_RE
      xbd%XT_XF = 0.0
C
      RETURN
      END

C===================================================================70
C 1st-order amplification equation
C===================================================================70
      SUBROUTINE TRCHEK(xfd,xbd)

      use xfoil_data_mod 
      type(xfoil_data_type), intent(inout) :: xfd
      type(xbl_data_type), intent(inout) :: xbd
C
cc      CALL TRCHEK1
C
C---- 2nd-order amplification equation
      CALL TRCHEK2(xfd,xbd)
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
      REAL MSQ
C
C---- density shape parameter    (from Whitfield)
      HC     = MSQ * (0.064/(HK-0.8) + 0.251)
      HC_HK  = MSQ * (-.064/(HK-0.8)**2     )
      HC_MSQ =        0.064/(HK-0.8) + 0.251
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      REAL MSQ
C
C---- Laminar HS correlation
      IF(HK.LT.4.35) THEN
       TMP = HK - 4.35
       HS    = 0.0111*TMP**2/(HK+1.0)
     &       - 0.0278*TMP**3/(HK+1.0)  + 1.528
     &       - 0.0002*(TMP*HK)**2
       HS_HK = 0.0111*(2.0*TMP    - TMP**2/(HK+1.0))/(HK+1.0)
     &       - 0.0278*(3.0*TMP**2 - TMP**3/(HK+1.0))/(HK+1.0)
     &       - 0.0002*2.0*TMP*HK * (TMP + HK)
      ELSE
       HS    = 0.015*    (HK-4.35)**2/HK + 1.528
       HS_HK = 0.015*2.0*(HK-4.35)   /HK
     &       - 0.015*    (HK-4.35)**2/HK**2
      ENDIF
C
      HS_RT  = 0.
      HS_MSQ = 0.
C
      RETURN
      END

C===================================================================70
C---- Turbulent HS correlation
C===================================================================70
      SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      IMPLICIT REAL (A-H,M,O-Z)
C
C
      DATA HSMIN, DHSINF / 1.500, 0.015 /
C
C---- ###  12/4/94
C---- limited Rtheta dependence for Rtheta < 200
C
C
      IF(RT.GT.400.0) THEN
       HO    = 3.0 + 400.0/RT
       HO_RT =     - 400.0/RT**2
      ELSE
       HO    = 4.0
       HO_RT = 0.
      ENDIF
C
      IF(RT.GT.200.0) THEN
       RTZ    = RT
       RTZ_RT = 1.
      ELSE
       RTZ    = 200.0
       RTZ_RT = 0.
      ENDIF
C
      IF(HK.LT.HO) THEN
C----- attached branch
C=======================================================
C----- old correlation
C-     (from Swafford profiles)
c       SRT = SQRT(RT)
c       HEX = (HO-HK)**1.6
c       RTMP = 0.165 - 1.6/SRT
c       HS    = HSMIN + 4.0/RT + RTMP*HEX/HK
c       HS_HK = RTMP*HEX/HK*(-1.6/(HO-HK) - 1.0/HK)
c       HS_RT = -4.0/RT**2 + HEX/HK*0.8/SRT/RT
c     &             + RTMP*HEX/HK*1.6/(HO-HK)*HO_RT
C=======================================================
C----- new correlation  29 Nov 91
C-     (from  arctan(y+) + Schlichting  profiles)
       HR    = ( HO - HK)/(HO-1.0)
       HR_HK =      - 1.0/(HO-1.0)
       HR_RT = (1.0 - HR)/(HO-1.0) * HO_RT
       HS    = (2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5) + HSMIN
     &       + 4.0/RTZ
       HS_HK =-(2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5)**2
     &       + (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_HK
       HS_RT = (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_RT
     &       + (HR**2 * 1.5/(HK+0.5) - 1.0)*4.0/RTZ**2 * RTZ_RT
C
      ELSE
C
C----- separated branch
       GRT = LOG(RTZ)
       HDIF = HK - HO 
       RTMP = HK - HO + 4.0/GRT
       HTMP    = 0.007*GRT/RTMP**2 + DHSINF/HK
       HTMP_HK = -.014*GRT/RTMP**3 - DHSINF/HK**2
       HTMP_RT = -.014*GRT/RTMP**3 * (-HO_RT - 4.0/GRT**2/RTZ * RTZ_RT)
     &         + 0.007    /RTMP**2 / RTZ * RTZ_RT
       HS    = HDIF**2 * HTMP + HSMIN + 4.0/RTZ
       HS_HK = HDIF*2.0* HTMP
     &       + HDIF**2 * HTMP_HK
       HS_RT = HDIF**2 * HTMP_RT      - 4.0/RTZ**2 * RTZ_RT
     &       + HDIF*2.0* HTMP * (-HO_RT)
C
      ENDIF
C
C---- fudge HS slightly to make sure   HS -> 2   as   HK -> 1
C-    (unnecessary with new correlation)
c      HTF    = 0.485/9.0 * (HK-4.0)**2/HK  +  1.515
c      HTF_HK = 0.485/9.0 * (1.0-16.0/HK**2)
c      ARG = MAX( 10.0*(1.0 - HK) , -15.0 )
c      HXX = EXP(ARG)
c      HXX_HK = -10.0*HXX
cC
c      HS_HK  = (1.0-HXX)*HS_HK  +  HXX*HTF_HK
c     &       + (        -HS     +      HTF    )*HXX_HK
c      HS_RT  = (1.0-HXX)*HS_RT
c      HS     = (1.0-HXX)*HS     +  HXX*HTF
C
C---- Whitfield's minor additional compressibility correction
      FM = 1.0 + 0.014*MSQ
      HS     = ( HS + 0.028*MSQ ) / FM
      HS_HK  = ( HS_HK          ) / FM
      HS_RT  = ( HS_RT          ) / FM
      HS_MSQ = 0.028/FM  -  0.014*HS/FM
C
      RETURN
      END

C===================================================================70
C---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
C===================================================================70
      SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      REAL*8 MSQ
C
      IF(HK.LT.5.5) THEN
       TMP = (5.5-HK)**3 / (HK+1.0)
       CF    = ( 0.0727*TMP                      - 0.07       )/RT
       CF_HK = ( -.0727*TMP*3.0/(5.5-HK) - 0.0727*TMP/(HK+1.0))/RT
      ELSE
       TMP = 1.0 - 1.0/(HK-4.5)
       CF    = ( 0.015*TMP**2      - 0.07  ) / RT
       CF_HK = ( 0.015*TMP*2.0/(HK-4.5)**2 ) / RT
      ENDIF
      CF_RT = -CF/RT
      CF_MSQ = 0.0
C
      RETURN
      END

C===================================================================70
C---- Turbulent skin friction function  ( Cf )    (Coles)
C===================================================================70
      SUBROUTINE CFT( bld, HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )

      use xfoil_data_mod
      IMPLICIT REAL (A-H,M,O-Z)
      type(blpar_data_type), intent(inout) :: bld
C
      DATA GAM /1.4/
C
      GM1 = GAM - 1.0
      FC = SQRT(1.0 + 0.5*GM1*MSQ)
      GRT = LOG(RT/FC)
      GRT = MAX(GRT,3.0)
C
      GEX = -1.74 - 0.31*HK
C
      ARG = -1.33*HK
      ARG = MAX(-20.0, ARG )
C
      THK = TANH(4.0 - HK/0.875)
C
      CFO =  bld%CFFAC * 0.3*EXP(ARG) * (GRT/2.3026)**GEX
      CF     = ( CFO  +  1.1E-4*(THK-1.0) ) / FC
      CF_HK  = (-1.33*CFO - 0.31*LOG(GRT/2.3026)*CFO
     &         - 1.1E-4*(1.0-THK**2) / 0.875    ) / FC
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.25*GM1/FC**2) - 0.25*GM1*CF/FC**2
C
      RETURN
      END ! CFT

C===================================================================70
C---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
C===================================================================70
      SUBROUTINE DIL( HK, RT, DI, DI_HK, DI_RT )
C
      IF(HK.LT.4.0) THEN
       DI    = ( 0.00205  *  (4.0-HK)**5.5 + 0.207 ) / RT
       DI_HK = ( -.00205*5.5*(4.0-HK)**4.5         ) / RT
      ELSE
       HKB = HK - 4.0
       DEN = 1.0 + 0.02*HKB**2
       DI    = ( -.0016  *  HKB**2  /DEN   + 0.207             ) / RT
       DI_HK = ( -.0016*2.0*HKB*(1.0/DEN - 0.02*HKB**2/DEN**2) ) / RT
      ENDIF
      DI_RT = -DI/RT
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE DILW( HK, RT, DI, DI_HK, DI_RT )
      REAL MSQ
C
      MSQ = 0.
      CALL HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
C
C---- Laminar wake dissipation function  ( 2 CD/H* )
      RCD    =  1.10 * (1.0 - 1.0/HK)**2  / HK
      RCD_HK = -1.10 * (1.0 - 1.0/HK)*2.0 / HK**3
     &       - RCD/HK
C
      DI    = 2.0*RCD   /(HS*RT)
      DI_HK = 2.0*RCD_HK/(HS*RT) - (DI/HS)*HS_HK
      DI_RT = -DI/RT             - (DI/HS)*HS_RT
C
      RETURN
      END

C===================================================================70
C
C     Calculates all secondary "2" variables from
C     the primary "2" variables X2, U2, T2, D2, S2.
C     Also calculates the sensitivities of the
C     secondary variables wrt the primary variables.
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C
C===================================================================70
      SUBROUTINE BLVAR(bld,xbd,ITYP)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
      IF(ITYP.EQ.3) xbd%HK2 = MAX(xbd%HK2,1.00005)
      IF(ITYP.NE.3) xbd%HK2 = MAX(xbd%HK2,1.05000)
C
C---- density thickness shape parameter     ( H** )
      CALL HCT( xbd%HK2, xbd%M2, xbd%HC2, HC2_HK2, HC2_M2 )
      xbd%HC2_U2 = HC2_HK2*xbd%HK2_U2 + HC2_M2*xbd%M2_U2
      xbd%HC2_T2 = HC2_HK2*xbd%HK2_T2
      xbd%HC2_D2 = HC2_HK2*xbd%HK2_D2
      xbd%HC2_MS = HC2_HK2*xbd%HK2_MS + HC2_M2*xbd%M2_MS
C
C---- set KE thickness shape parameter from  H - H*  correlations
      IF(ITYP.EQ.1) THEN
       CALL HSL( xbd%HK2, xbd%RT2, xbd%M2, xbd%HS2, HS2_HK2, HS2_RT2,
     &   HS2_M2 )
      ELSE
       CALL HST( xbd%HK2, xbd%RT2, xbd%M2, xbd%HS2, HS2_HK2, HS2_RT2,
     &   HS2_M2 )
      ENDIF
C
      xbd%HS2_U2 = HS2_HK2*xbd%HK2_U2 + HS2_RT2*xbd%RT2_U2 + HS2_M2
     &  *xbd%M2_U2
      xbd%HS2_T2 = HS2_HK2*xbd%HK2_T2 + HS2_RT2*xbd%RT2_T2
      xbd%HS2_D2 = HS2_HK2*xbd%HK2_D2
      xbd%HS2_MS = HS2_HK2*xbd%HK2_MS + HS2_RT2*xbd%RT2_MS + HS2_M2
     &  *xbd%M2_MS
      xbd%HS2_RE =                  HS2_RT2*xbd%RT2_RE
C
C---- normalized slip velocity  Us
      xbd%US2     = 0.5*xbd%HS2*( 1.0 - (xbd%HK2-1.0)/(bld%GBCON*xbd%H2)
     &   )
      US2_HS2 = 0.5  *  ( 1.0 - (xbd%HK2-1.0)/(bld%GBCON*xbd%H2) )
      US2_HK2 = 0.5*xbd%HS2*(     -  1.0     /(bld%GBCON*xbd%H2) )
      US2_H2  = 0.5*xbd%HS2*        (xbd%HK2-1.0)/(bld%GBCON*xbd%H2**2)
C
      xbd%US2_U2 = US2_HS2*xbd%HS2_U2 + US2_HK2*xbd%HK2_U2
      xbd%US2_T2 = US2_HS2*xbd%HS2_T2 + US2_HK2*xbd%HK2_T2 + US2_H2
     &  *xbd%H2_T2
      xbd%US2_D2 = US2_HS2*xbd%HS2_D2 + US2_HK2*xbd%HK2_D2 + US2_H2
     &  *xbd%H2_D2
      xbd%US2_MS = US2_HS2*xbd%HS2_MS + US2_HK2*xbd%HK2_MS
      xbd%US2_RE = US2_HS2*xbd%HS2_RE
C
      IF(ITYP.LE.2 .AND. xbd%US2.GT.0.95) THEN
CCC       WRITE(*,*) 'BLVAR: Us clamped:', US2
       xbd%US2 = 0.98
       xbd%US2_U2 = 0.
       xbd%US2_T2 = 0.
       xbd%US2_D2 = 0.
       xbd%US2_MS = 0.
       xbd%US2_RE = 0.
      ENDIF
C
      IF(ITYP.EQ.3 .AND. xbd%US2.GT.0.99995) THEN
CCC       WRITE(*,*) 'BLVAR: Wake Us clamped:', US2
       xbd%US2 = 0.99995
       xbd%US2_U2 = 0.
       xbd%US2_T2 = 0.
       xbd%US2_D2 = 0.
       xbd%US2_MS = 0.
       xbd%US2_RE = 0.
      ENDIF
C
C---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
C   ...  NEW  12 Oct 94
      GCC = 0.0
      HKC = xbd%HK2 - 1.0
      HKC_HK2 = 1.0
      HKC_RT2 = 0.0
      IF(ITYP.EQ.2) THEN
       GCC = bld%GCCON
       HKC     = xbd%HK2 - 1.0 - GCC/xbd%RT2
       HKC_HK2 = 1.0
       HKC_RT2 =             GCC/xbd%RT2**2
       IF(HKC .LT. 0.01) THEN
        HKC = 0.01
        HKC_HK2 = 0.0
        HKC_RT2 = 0.0
       ENDIF
      ENDIF
C
      HKB = xbd%HK2 - 1.0
      USB = 1.0 - xbd%US2
      xbd%CQ2     =
     &    SQRT( bld%CTCON*xbd%HS2*HKB*HKC**2 / (USB*xbd%H2*xbd%HK2**2) )
     &  
      CQ2_HS2 = bld%CTCON    *HKB*HKC**2 / (USB*xbd%H2*xbd%HK2**2)      
     &   * 0.5/xbd%CQ2
      CQ2_US2 = bld%CTCON*xbd%HS2*HKB*HKC**2 / (USB*xbd%H2*xbd%HK2**2) /
     &   USB * 0.5/xbd%CQ2
      CQ2_HK2 = bld%CTCON*xbd%HS2    *HKC**2 / (USB*xbd%H2*xbd%HK2**2)  
     &       * 0.5/xbd%CQ2
     &        - bld%CTCON*xbd%HS2*HKB*HKC**2 / (USB*xbd%H2*xbd%HK2**3) *
     &   2.0 * 0.5/xbd%CQ2
     &        + bld%CTCON*xbd%HS2*HKB*HKC    / (USB*xbd%H2*xbd%HK2**2) *
     &   2.0 * 0.5/xbd%CQ2
     &         *HKC_HK2
      CQ2_RT2 = bld%CTCON*xbd%HS2*HKB*HKC    / (USB*xbd%H2*xbd%HK2**2) *
     &   2.0 * 0.5/xbd%CQ2
     &         *HKC_RT2
      CQ2_H2  =-bld%CTCON*xbd%HS2*HKB*HKC**2 / (USB*xbd%H2*xbd%HK2**2) /
     &   xbd%H2  * 0.5/xbd%CQ2
C
      xbd%CQ2_U2 = CQ2_HS2*xbd%HS2_U2 + CQ2_US2*xbd%US2_U2 + CQ2_HK2
     &  *xbd%HK2_U2
      xbd%CQ2_T2 = CQ2_HS2*xbd%HS2_T2 + CQ2_US2*xbd%US2_T2 + CQ2_HK2
     &  *xbd%HK2_T2
      xbd%CQ2_D2 = CQ2_HS2*xbd%HS2_D2 + CQ2_US2*xbd%US2_D2 + CQ2_HK2
     &  *xbd%HK2_D2
      xbd%CQ2_MS = CQ2_HS2*xbd%HS2_MS + CQ2_US2*xbd%US2_MS + CQ2_HK2
     &  *xbd%HK2_MS
      xbd%CQ2_RE = CQ2_HS2*xbd%HS2_RE + CQ2_US2*xbd%US2_RE
C
      xbd%CQ2_U2 = xbd%CQ2_U2                + CQ2_RT2*xbd%RT2_U2
      xbd%CQ2_T2 = xbd%CQ2_T2 + CQ2_H2*xbd%H2_T2 + CQ2_RT2*xbd%RT2_T2
      xbd%CQ2_D2 = xbd%CQ2_D2 + CQ2_H2*xbd%H2_D2
      xbd%CQ2_MS = xbd%CQ2_MS                + CQ2_RT2*xbd%RT2_MS
      xbd%CQ2_RE = xbd%CQ2_RE                + CQ2_RT2*xbd%RT2_RE
C
C
C---- set skin friction coefficient 
      IF(ITYP.EQ.3) THEN
C----- wake
       xbd%CF2     = 0.
       CF2_HK2 = 0.
       CF2_RT2 = 0.
       CF2_M2  = 0.
      ELSE IF(ITYP.EQ.1) THEN
C----- laminar
       CALL CFL( xbd%HK2, xbd%RT2, xbd%M2, xbd%CF2, CF2_HK2, CF2_RT2,
     &   CF2_M2 )
      ELSE
C----- turbulent
       CALL CFT( bld, xbd%HK2, xbd%RT2, xbd%M2, xbd%CF2, CF2_HK2,
     &   CF2_RT2, CF2_M2 )
       CALL CFL( xbd%HK2, xbd%RT2, xbd%M2, CF2L,CF2L_HK2,CF2L_RT2
     &  ,CF2L_M2)
       IF(CF2L.GT.xbd%CF2) THEN
C------- laminar Cf is greater than turbulent Cf -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc      write(*,*) 'Cft Cfl Rt Hk:', CF2, CF2L, RT2, HK2, X2
         xbd%CF2     = CF2L
         CF2_HK2 = CF2L_HK2
         CF2_RT2 = CF2L_RT2
         CF2_M2  = CF2L_M2
       ENDIF
      ENDIF
C
      xbd%CF2_U2 = CF2_HK2*xbd%HK2_U2 + CF2_RT2*xbd%RT2_U2 + CF2_M2
     &  *xbd%M2_U2
      xbd%CF2_T2 = CF2_HK2*xbd%HK2_T2 + CF2_RT2*xbd%RT2_T2
      xbd%CF2_D2 = CF2_HK2*xbd%HK2_D2
      xbd%CF2_MS = CF2_HK2*xbd%HK2_MS + CF2_RT2*xbd%RT2_MS + CF2_M2
     &  *xbd%M2_MS
      xbd%CF2_RE =                  CF2_RT2*xbd%RT2_RE
C
C---- dissipation function    2 CD / H*
      IF(ITYP.EQ.1) THEN
C
C----- laminar
       CALL DIL( xbd%HK2, xbd%RT2, xbd%DI2, DI2_HK2, DI2_RT2 )
C
       xbd%DI2_U2 = DI2_HK2*xbd%HK2_U2 + DI2_RT2*xbd%RT2_U2
       xbd%DI2_T2 = DI2_HK2*xbd%HK2_T2 + DI2_RT2*xbd%RT2_T2
       xbd%DI2_D2 = DI2_HK2*xbd%HK2_D2
       xbd%DI2_S2 = 0.
       xbd%DI2_MS = DI2_HK2*xbd%HK2_MS + DI2_RT2*xbd%RT2_MS
       xbd%DI2_RE =                  DI2_RT2*xbd%RT2_RE
C
      ELSE IF(ITYP.EQ.2) THEN
C
CCC       CALL DIT(     HS2,     US2,     CF2,     S2, DI2,
CCC     &           DI2_HS2, DI2_US2, DI2_CF2, DI2_S2      )
C
C----- turbulent wall contribution
       CALL CFT(bld, xbd%HK2, xbd%RT2, xbd%M2, CF2T, CF2T_HK2, CF2T_RT2,
     &   CF2T_M2)
       CF2T_U2 = CF2T_HK2*xbd%HK2_U2 + CF2T_RT2*xbd%RT2_U2 + CF2T_M2
     &  *xbd%M2_U2
       CF2T_T2 = CF2T_HK2*xbd%HK2_T2 + CF2T_RT2*xbd%RT2_T2
       CF2T_D2 = CF2T_HK2*xbd%HK2_D2
       CF2T_MS = CF2T_HK2*xbd%HK2_MS + CF2T_RT2*xbd%RT2_MS + CF2T_M2
     &  *xbd%M2_MS
       CF2T_RE =                   CF2T_RT2*xbd%RT2_RE
C
       xbd%DI2      =  ( 0.5*CF2T*xbd%US2 ) * 2.0/xbd%HS2
       DI2_HS2  = -( 0.5*CF2T*xbd%US2 ) * 2.0/xbd%HS2**2
       DI2_US2  =  ( 0.5*CF2T     ) * 2.0/xbd%HS2
       DI2_CF2T =  ( 0.5     *xbd%US2 ) * 2.0/xbd%HS2
C
       xbd%DI2_S2 = 0.0
       xbd%DI2_U2 = DI2_HS2*xbd%HS2_U2 + DI2_US2*xbd%US2_U2 + DI2_CF2T
     &  *CF2T_U2
       xbd%DI2_T2 = DI2_HS2*xbd%HS2_T2 + DI2_US2*xbd%US2_T2 + DI2_CF2T
     &  *CF2T_T2
       xbd%DI2_D2 = DI2_HS2*xbd%HS2_D2 + DI2_US2*xbd%US2_D2 + DI2_CF2T
     &  *CF2T_D2
       xbd%DI2_MS = DI2_HS2*xbd%HS2_MS + DI2_US2*xbd%US2_MS + DI2_CF2T
     &  *CF2T_MS
       xbd%DI2_RE = DI2_HS2*xbd%HS2_RE + DI2_US2*xbd%US2_RE + DI2_CF2T
     &  *CF2T_RE
C
C
C----- set minimum Hk for wake layer to still exist
       GRT = LOG(xbd%RT2)
       HMIN = 1.0 + 2.1/GRT
       HM_RT2 = -(2.1/GRT**2) / xbd%RT2
C
C----- set factor DFAC for correcting wall dissipation for very low Hk
       FL = (xbd%HK2-1.0)/(HMIN-1.0)
       FL_HK2 =   1.0/(HMIN-1.0)
       FL_RT2 = ( -FL/(HMIN-1.0) ) * HM_RT2
C
       TFL = TANH(FL)
       DFAC  = 0.5 + 0.5* TFL
       DF_FL =       0.5*(1.0 - TFL**2)
C
       DF_HK2 = DF_FL*FL_HK2
       DF_RT2 = DF_FL*FL_RT2
C
       xbd%DI2_S2 = xbd%DI2_S2*DFAC
       xbd%DI2_U2 = xbd%DI2_U2*DFAC + xbd%DI2*(DF_HK2*xbd%HK2_U2 +
     &   DF_RT2*xbd%RT2_U2)
       xbd%DI2_T2 = xbd%DI2_T2*DFAC + xbd%DI2*(DF_HK2*xbd%HK2_T2 +
     &   DF_RT2*xbd%RT2_T2)
       xbd%DI2_D2 = xbd%DI2_D2*DFAC + xbd%DI2*(DF_HK2*xbd%HK2_D2        
     &          )
       xbd%DI2_MS = xbd%DI2_MS*DFAC + xbd%DI2*(DF_HK2*xbd%HK2_MS +
     &   DF_RT2*xbd%RT2_MS)
       xbd%DI2_RE = xbd%DI2_RE*DFAC + xbd%DI2*(                DF_RT2
     &  *xbd%RT2_RE)
       xbd%DI2    = xbd%DI2   *DFAC
C
      ELSE
C
C----- zero wall contribution for wake
       xbd%DI2    = 0.0
       xbd%DI2_S2 = 0.0
       xbd%DI2_U2 = 0.0
       xbd%DI2_T2 = 0.0
       xbd%DI2_D2 = 0.0
       xbd%DI2_MS = 0.0
       xbd%DI2_RE = 0.0
C
      ENDIF
C
C
C---- Add on turbulent outer layer contribution
      IF(ITYP.NE.1) THEN
C
       DD     =  xbd%S2**2 * (0.995-xbd%US2) * 2.0/xbd%HS2
       DD_HS2 = -xbd%S2**2 * (0.995-xbd%US2) * 2.0/xbd%HS2**2
       DD_US2 = -xbd%S2**2               * 2.0/xbd%HS2
       DD_S2  =  xbd%S2*2.0* (0.995-xbd%US2) * 2.0/xbd%HS2
C
       xbd%DI2    = xbd%DI2    + DD
       xbd%DI2_S2 =          DD_S2
       xbd%DI2_U2 = xbd%DI2_U2 + DD_HS2*xbd%HS2_U2 + DD_US2*xbd%US2_U2
       xbd%DI2_T2 = xbd%DI2_T2 + DD_HS2*xbd%HS2_T2 + DD_US2*xbd%US2_T2
       xbd%DI2_D2 = xbd%DI2_D2 + DD_HS2*xbd%HS2_D2 + DD_US2*xbd%US2_D2
       xbd%DI2_MS = xbd%DI2_MS + DD_HS2*xbd%HS2_MS + DD_US2*xbd%US2_MS
       xbd%DI2_RE = xbd%DI2_RE + DD_HS2*xbd%HS2_RE + DD_US2*xbd%US2_RE
C
C----- add laminar stress contribution to outer layer CD
c###
       DD     =  0.15*(0.995-xbd%US2)**2 / xbd%RT2  * 2.0/xbd%HS2
       DD_US2 = -0.15*(0.995-xbd%US2)*2. / xbd%RT2  * 2.0/xbd%HS2
       DD_HS2 = -DD/xbd%HS2
       DD_RT2 = -DD/xbd%RT2
C
       xbd%DI2    = xbd%DI2    + DD
       xbd%DI2_U2 = xbd%DI2_U2 + DD_HS2*xbd%HS2_U2 + DD_US2*xbd%US2_U2 +
     &   DD_RT2*xbd%RT2_U2
       xbd%DI2_T2 = xbd%DI2_T2 + DD_HS2*xbd%HS2_T2 + DD_US2*xbd%US2_T2 +
     &   DD_RT2*xbd%RT2_T2
       xbd%DI2_D2 = xbd%DI2_D2 + DD_HS2*xbd%HS2_D2 + DD_US2*xbd%US2_D2
       xbd%DI2_MS = xbd%DI2_MS + DD_HS2*xbd%HS2_MS + DD_US2*xbd%US2_MS +
     &   DD_RT2*xbd%RT2_MS
       xbd%DI2_RE = xbd%DI2_RE + DD_HS2*xbd%HS2_RE + DD_US2*xbd%US2_RE +
     &   DD_RT2*xbd%RT2_RE
C
      ENDIF
C
C
      IF(ITYP.EQ.2) THEN
        CALL DIL( xbd%HK2, xbd%RT2, DI2L, DI2L_HK2, DI2L_RT2 )
C
        IF(DI2L.GT.xbd%DI2) THEN
C------- laminar CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc       write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
          xbd%DI2    = DI2L
          xbd%DI2_S2 = 0.
          xbd%DI2_U2 = DI2L_HK2*xbd%HK2_U2 + DI2L_RT2*xbd%RT2_U2
          xbd%DI2_T2 = DI2L_HK2*xbd%HK2_T2 + DI2L_RT2*xbd%RT2_T2
          xbd%DI2_D2 = DI2L_HK2*xbd%HK2_D2
          xbd%DI2_MS = DI2L_HK2*xbd%HK2_MS + DI2L_RT2*xbd%RT2_MS
          xbd%DI2_RE =                   DI2L_RT2*xbd%RT2_RE
        ENDIF
      ENDIF
C
cC----- add on CD contribution of inner shear layer
c       IF(ITYP.EQ.3 .AND. DW2.GT.0.0) THEN
c        DKON = 0.03*0.75**3
c        DDI = DKON*US2**3
c        DDI_US2 = 3.0*DKON*US2**2
c        DI2 = DI2 + DDI * DW2/DWTE
c        DI2_U2 = DI2_U2 + DDI_US2*US2_U2 * DW2/DWTE
c        DI2_T2 = DI2_T2 + DDI_US2*US2_T2 * DW2/DWTE
c        DI2_D2 = DI2_D2 + DDI_US2*US2_D2 * DW2/DWTE
c        DI2_MS = DI2_MS + DDI_US2*US2_MS * DW2/DWTE
c        DI2_RE = DI2_RE + DDI_US2*US2_RE * DW2/DWTE
c       ENDIF
C
      IF(ITYP.EQ.3) THEN
C------ laminar wake CD
        CALL DILW( xbd%HK2, xbd%RT2, DI2L, DI2L_HK2, DI2L_RT2 )
        IF(DI2L .GT. xbd%DI2) THEN
C------- laminar wake CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc         write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
         xbd%DI2    = DI2L
         xbd%DI2_S2 = 0.
         xbd%DI2_U2 = DI2L_HK2*xbd%HK2_U2 + DI2L_RT2*xbd%RT2_U2
         xbd%DI2_T2 = DI2L_HK2*xbd%HK2_T2 + DI2L_RT2*xbd%RT2_T2
         xbd%DI2_D2 = DI2L_HK2*xbd%HK2_D2
         xbd%DI2_MS = DI2L_HK2*xbd%HK2_MS + DI2L_RT2*xbd%RT2_MS
         xbd%DI2_RE =                   DI2L_RT2*xbd%RT2_RE
        ENDIF
      ENDIF
C
C
      IF(ITYP.EQ.3) THEN
C----- double dissipation for the wake (two wake halves)
       xbd%DI2    = xbd%DI2   *2.0
       xbd%DI2_S2 = xbd%DI2_S2*2.0
       xbd%DI2_U2 = xbd%DI2_U2*2.0
       xbd%DI2_T2 = xbd%DI2_T2*2.0
       xbd%DI2_D2 = xbd%DI2_D2*2.0
       xbd%DI2_MS = xbd%DI2_MS*2.0
       xbd%DI2_RE = xbd%DI2_RE*2.0
      ENDIF
C
C---- BL thickness (Delta) from simplified Green's correlation
      xbd%DE2     = (3.15 + 1.72/(xbd%HK2-1.0)   )*xbd%T2  +  xbd%D2
      DE2_HK2 = (     - 1.72/(xbd%HK2-1.0)**2)*xbd%T2
C
      xbd%DE2_U2 = DE2_HK2*xbd%HK2_U2
      xbd%DE2_T2 = DE2_HK2*xbd%HK2_T2 + (3.15 + 1.72/(xbd%HK2-1.0))
      xbd%DE2_D2 = DE2_HK2*xbd%HK2_D2 + 1.0
      xbd%DE2_MS = DE2_HK2*xbd%HK2_MS
C
ccc      HDMAX = 15.0
      HDMAX = 12.0
      IF(xbd%DE2 .GT. HDMAX*xbd%T2) THEN
cccc      IF(DE2 .GT. HDMAX*T2 .AND. (HK2 .GT. 4.0 .OR. ITYP.EQ.3)) THEN
       xbd%DE2    = HDMAX*xbd%T2
       xbd%DE2_U2 =  0.0
       xbd%DE2_T2 = HDMAX
       xbd%DE2_D2 =  0.0
       xbd%DE2_MS =  0.0
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Sets up "dummy" BL system between airfoil TE point 
C     and first wake point infinitesimally behind TE.
C
C===================================================================70
      SUBROUTINE TESYS(bld,xbd,CTE,TTE,DTE)

      use xfoil_data_mod
      IMPLICIT REAL (M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
      DO 55 K=1, 4
        xbd%VSREZ(K) = 0.
        xbd%VSM(K)   = 0.
        xbd%VSR(K)   = 0.
        xbd%VSX(K)   = 0.
        DO 551 L=1, 5
          xbd%VS1(K,L) = 0.
          xbd%VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
      CALL BLVAR(bld,xbd,3)
C
      xbd%VS1(1,1) = -1.0
      xbd%VS2(1,1) = 1.0
      xbd%VSREZ(1) = CTE - xbd%S2      
C
      xbd%VS1(2,2) = -1.0
      xbd%VS2(2,2) = 1.0
      xbd%VSREZ(2) = TTE - xbd%T2
C
      xbd%VS1(3,3) = -1.0
      xbd%VS2(3,3) = 1.0
      xbd%VSREZ(3) = DTE - xbd%D2 - xbd%DW2
C
      RETURN
      END

C===================================================================70
C
C     Calculates midpoint skin friction CFM
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C
C===================================================================70
      SUBROUTINE BLMID(bld,xbd,ITYP)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
C---- set similarity variables if not defined
      IF(xbd%SIMI) THEN
       xbd%HK1    = xbd%HK2
       xbd%HK1_T1 = xbd%HK2_T2
       xbd%HK1_D1 = xbd%HK2_D2
       xbd%HK1_U1 = xbd%HK2_U2
       xbd%HK1_MS = xbd%HK2_MS
       xbd%RT1    = xbd%RT2
       xbd%RT1_T1 = xbd%RT2_T2
       xbd%RT1_U1 = xbd%RT2_U2
       xbd%RT1_MS = xbd%RT2_MS
       xbd%RT1_RE = xbd%RT2_RE
       xbd%M1    = xbd%M2
       xbd%M1_U1 = xbd%M2_U2
       xbd%M1_MS = xbd%M2_MS
      ENDIF
C
C---- define stuff for midpoint CF
      HKA = 0.5*(xbd%HK1 + xbd%HK2)
      RTA = 0.5*(xbd%RT1 + xbd%RT2)
      MA  = 0.5*(xbd%M1  + xbd%M2 )
C
C---- midpoint skin friction coefficient  (zero in wake)
      IF(ITYP.EQ.3) THEN
       xbd%CFM     = 0.
       CFM_HKA = 0.
       CFM_RTA = 0.
       CFM_MA  = 0.
       xbd%CFM_MS  = 0.
      ELSE IF(ITYP.EQ.1) THEN
       CALL CFL( HKA, RTA, MA, xbd%CFM, CFM_HKA, CFM_RTA, CFM_MA )
      ELSE
       CALL CFT( bld, HKA, RTA, MA, xbd%CFM, CFM_HKA, CFM_RTA, CFM_MA )
       CALL CFL( HKA, RTA, MA, CFML,CFML_HKA,CFML_RTA,CFML_MA)
       IF(CFML.GT.xbd%CFM) THEN
ccc      write(*,*) 'Cft Cfl Rt Hk:', CFM, CFML, RTA, HKA, 0.5*(X1+X2)
         xbd%CFM     = CFML
         CFM_HKA = CFML_HKA
         CFM_RTA = CFML_RTA
         CFM_MA  = CFML_MA
       ENDIF
      ENDIF
C
      xbd%CFM_U1 = 0.5*(CFM_HKA*xbd%HK1_U1 + CFM_MA*xbd%M1_U1 + CFM_RTA
     &  *xbd%RT1_U1)
      xbd%CFM_T1 = 0.5*(CFM_HKA*xbd%HK1_T1 +                CFM_RTA
     &  *xbd%RT1_T1)
      xbd%CFM_D1 = 0.5*(CFM_HKA*xbd%HK1_D1                              
     &    )
C
      xbd%CFM_U2 = 0.5*(CFM_HKA*xbd%HK2_U2 + CFM_MA*xbd%M2_U2 + CFM_RTA
     &  *xbd%RT2_U2)
      xbd%CFM_T2 = 0.5*(CFM_HKA*xbd%HK2_T2 +                CFM_RTA
     &  *xbd%RT2_T2)
      xbd%CFM_D2 = 0.5*(CFM_HKA*xbd%HK2_D2                              
     &    )
C
      xbd%CFM_MS = 0.5*(CFM_HKA*xbd%HK1_MS + CFM_MA*xbd%M1_MS + CFM_RTA
     &  *xbd%RT1_MS
     &            + CFM_HKA*xbd%HK2_MS + CFM_MA*xbd%M2_MS + CFM_RTA
     &  *xbd%RT2_MS)
      xbd%CFM_RE = 0.5*(                                CFM_RTA
     &  *xbd%RT1_RE
     &                                            + CFM_RTA*xbd%RT2_RE)
C
      RETURN
      END ! BLMID

C===================================================================70
C
C     Sets up the Newton system coefficients and residuals
C
C        ITYP = 0 :  similarity station
C        ITYP = 1 :  laminar interval
C        ITYP = 2 :  turbulent interval
C        ITYP = 3 :  wake interval
C
C     This routine knows nothing about a transition interval,
C     which is taken care of by TRDIF.
C
C===================================================================70
      SUBROUTINE BLDIF(bld,xbd,ITYP)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
      IF(ITYP.EQ.0) THEN
C----- similarity logarithmic differences  (prescribed)
       XLOG = 1.0
       ULOG = xbd%BULE
       TLOG = 0.5*(1.0 - xbd%BULE)
       HLOG = 0.
       DDLOG = 0.
      ELSE
C----- usual logarithmic differences
       XLOG = LOG(xbd%X2/xbd%X1)
       ULOG = LOG(xbd%U2/xbd%U1)
       TLOG = LOG(xbd%T2/xbd%T1)
       HLOG = LOG(xbd%HS2/xbd%HS1)
C       XLOG = 2.0*(X2-X1)/(X2+X1)
C       ULOG = 2.0*(U2-U1)/(U2+U1)
C       TLOG = 2.0*(T2-T1)/(T2+T1)
C       HLOG = 2.0*(HS2-HS1)/(HS2+HS1)
       DDLOG = 1.0
      ENDIF
C
      DO 55 K=1, 4
        xbd%VSREZ(K) = 0.
        xbd%VSM(K) = 0.
        xbd%VSR(K) = 0.
        xbd%VSX(K) = 0.
        DO 551 L=1, 5
          xbd%VS1(K,L) = 0.
          xbd%VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
C---- set triggering constant for local upwinding
      HUPWT = 1.0
C
ccc      HDCON = 5.0*HUPWT
ccc      HD_HK1 = 0.0
ccc      HD_HK2 = 0.0
C
      HDCON  =  5.0*HUPWT/xbd%HK2**2
      HD_HK1 =  0.0
      HD_HK2 = -HDCON*2.0/xbd%HK2
C
C---- use less upwinding in the wake
      IF(ITYP.EQ.3) THEN
       HDCON  =  HUPWT/xbd%HK2**2
       HD_HK1 =  0.0
       HD_HK2 = -HDCON*2.0/xbd%HK2
      ENDIF
C
C---- local upwinding is based on local change in  log(Hk-1)
C-    (mainly kicks in at transition)
      ARG = ABS((xbd%HK2-1.0)/(xbd%HK1-1.0))
      HL = LOG(ARG)
      HL_HK1 = -1.0/(xbd%HK1-1.0)
      HL_HK2 =  1.0/(xbd%HK2-1.0)
C
C---- set local upwinding parameter UPW and linearize it
C
C       UPW = 0.5   Trapezoidal
C       UPW = 1.0   Backward Euler
C
      HLSQ = MIN( HL**2 , 15.0 )
      EHH = EXP(-HLSQ*HDCON)
      UPW = 1.0 - 0.5*EHH
      UPW_HL =        EHH * HL  *HDCON
      UPW_HD =    0.5*EHH * HLSQ
C
      UPW_HK1 = UPW_HL*HL_HK1 + UPW_HD*HD_HK1
      UPW_HK2 = UPW_HL*HL_HK2 + UPW_HD*HD_HK2
C
      UPW_U1 = UPW_HK1*xbd%HK1_U1
      UPW_T1 = UPW_HK1*xbd%HK1_T1
      UPW_D1 = UPW_HK1*xbd%HK1_D1
      UPW_U2 = UPW_HK2*xbd%HK2_U2
      UPW_T2 = UPW_HK2*xbd%HK2_T2
      UPW_D2 = UPW_HK2*xbd%HK2_D2
      UPW_MS = UPW_HK1*xbd%HK1_MS
     &       + UPW_HK2*xbd%HK2_MS
C
C
      IF(ITYP.EQ.0) THEN
C
C***** LE point -->  set zero amplification factor
       xbd%VS2(1,1) = 1.0
       xbd%VSR(1)   = 0.
       xbd%VSREZ(1) = -xbd%AMPL2
C
      ELSE IF(ITYP.EQ.1) THEN
C
C***** laminar part -->  set amplification equation
C
C----- set average amplification AX over interval X1..X2
       CALL AXSET( xbd%HK1,    xbd%T1,    xbd%RT1, xbd%AMPL1,  
     &             xbd%HK2,    xbd%T2,    xbd%RT2, xbd%AMPL2, xbd%AMCRIT
     &  , xbd%IDAMPV,
     &      AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &          AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
       REZC = xbd%AMPL2 - xbd%AMPL1 - AX*(xbd%X2-xbd%X1)
       Z_AX = -(xbd%X2-xbd%X1)
C
       xbd%VS1(1,1) = Z_AX* AX_A1  -  1.0
       xbd%VS1(1,2) = Z_AX*(AX_HK1*xbd%HK1_T1 + AX_T1 + AX_RT1
     &  *xbd%RT1_T1)
       xbd%VS1(1,3) = Z_AX*(AX_HK1*xbd%HK1_D1                        )
       xbd%VS1(1,4) = Z_AX*(AX_HK1*xbd%HK1_U1         + AX_RT1
     &  *xbd%RT1_U1)
       xbd%VS1(1,5) =  AX
       xbd%VS2(1,1) = Z_AX* AX_A2  +  1.0
       xbd%VS2(1,2) = Z_AX*(AX_HK2*xbd%HK2_T2 + AX_T2 + AX_RT2
     &  *xbd%RT2_T2)
       xbd%VS2(1,3) = Z_AX*(AX_HK2*xbd%HK2_D2                        )  
     &         
       xbd%VS2(1,4) = Z_AX*(AX_HK2*xbd%HK2_U2         + AX_RT2
     &  *xbd%RT2_U2)
       xbd%VS2(1,5) = -AX
       xbd%VSM(1)   = Z_AX*(AX_HK1*xbd%HK1_MS         + AX_RT1
     &  *xbd%RT1_MS
     &                + AX_HK2*xbd%HK2_MS         + AX_RT2*xbd%RT2_MS)
       xbd%VSR(1)   = Z_AX*(                        AX_RT1*xbd%RT1_RE
     &                                        + AX_RT2*xbd%RT2_RE)
       xbd%VSX(1)   = 0.
       xbd%VSREZ(1) = -REZC
C
      ELSE
C
C***** turbulent part -->  set shear lag equation
C
       SA  = (1.0-UPW)*xbd%S1  + UPW*xbd%S2
       CQA = (1.0-UPW)*xbd%CQ1 + UPW*xbd%CQ2
       CFA = (1.0-UPW)*xbd%CF1 + UPW*xbd%CF2
       HKA = (1.0-UPW)*xbd%HK1 + UPW*xbd%HK2
C
       USA = 0.5*(xbd%US1 + xbd%US2)
       RTA = 0.5*(xbd%RT1 + xbd%RT2)
       DEA = 0.5*(xbd%DE1 + xbd%DE2)
       DA  = 0.5*(xbd%D1  + xbd%D2 )
C
C
       IF(ITYP.EQ.3) THEN
C------ increased dissipation length in wake (decrease its reciprocal)
        ALD = bld%DLCON
       ELSE
        ALD = 1.0
       ENDIF
C
C----- set and linearize  equilibrium 1/Ue dUe/dx   ...  NEW  12 Oct 94
       IF(ITYP.EQ.2) THEN
        GCC = bld%GCCON
        HKC     = HKA - 1.0 - GCC/RTA
        HKC_HKA = 1.0
        HKC_RTA =             GCC/RTA**2
        IF(HKC .LT. 0.01) THEN
         HKC = 0.01
         HKC_HKA = 0.0
         HKC_RTA = 0.0
        ENDIF
       ELSE
        GCC = 0.0
        HKC = HKA - 1.0
        HKC_HKA = 1.0
        HKC_RTA = 0.0
       ENDIF
C
       HR     = HKC     / (bld%GACON*ALD*HKA)
       HR_HKA = HKC_HKA / (bld%GACON*ALD*HKA) - HR / HKA
       HR_RTA = HKC_RTA / (bld%GACON*ALD*HKA)
C
       UQ     = (0.5*CFA - HR**2) / (bld%GBCON*DA)
       UQ_HKA =   -2.0*HR*HR_HKA  / (bld%GBCON*DA)
       UQ_RTA =   -2.0*HR*HR_RTA  / (bld%GBCON*DA)
       UQ_CFA =   0.5             / (bld%GBCON*DA)
       UQ_DA  = -UQ/DA
       UQ_UPW = UQ_CFA*(xbd%CF2-xbd%CF1) + UQ_HKA*(xbd%HK2-xbd%HK1)
C
       UQ_T1 = (1.0-UPW)*(UQ_CFA*xbd%CF1_T1 + UQ_HKA*xbd%HK1_T1) +
     &   UQ_UPW*UPW_T1
       UQ_D1 = (1.0-UPW)*(UQ_CFA*xbd%CF1_D1 + UQ_HKA*xbd%HK1_D1) +
     &   UQ_UPW*UPW_D1
       UQ_U1 = (1.0-UPW)*(UQ_CFA*xbd%CF1_U1 + UQ_HKA*xbd%HK1_U1) +
     &   UQ_UPW*UPW_U1
       UQ_T2 =      UPW *(UQ_CFA*xbd%CF2_T2 + UQ_HKA*xbd%HK2_T2) +
     &   UQ_UPW*UPW_T2
       UQ_D2 =      UPW *(UQ_CFA*xbd%CF2_D2 + UQ_HKA*xbd%HK2_D2) +
     &   UQ_UPW*UPW_D2
       UQ_U2 =      UPW *(UQ_CFA*xbd%CF2_U2 + UQ_HKA*xbd%HK2_U2) +
     &   UQ_UPW*UPW_U2
       UQ_MS = (1.0-UPW)*(UQ_CFA*xbd%CF1_MS + UQ_HKA*xbd%HK1_MS) +
     &   UQ_UPW*UPW_MS
     &       +      UPW *(UQ_CFA*xbd%CF2_MS + UQ_HKA*xbd%HK2_MS)
       UQ_RE = (1.0-UPW)* UQ_CFA*xbd%CF1_RE
     &       +      UPW * UQ_CFA*xbd%CF2_RE
C
       UQ_T1 = UQ_T1             + 0.5*UQ_RTA*xbd%RT1_T1
       UQ_D1 = UQ_D1 + 0.5*UQ_DA
       UQ_U1 = UQ_U1             + 0.5*UQ_RTA*xbd%RT1_U1
       UQ_T2 = UQ_T2             + 0.5*UQ_RTA*xbd%RT2_T2
       UQ_D2 = UQ_D2 + 0.5*UQ_DA
       UQ_U2 = UQ_U2             + 0.5*UQ_RTA*xbd%RT2_U2
       UQ_MS = UQ_MS             + 0.5*UQ_RTA*xbd%RT1_MS
     &                           + 0.5*UQ_RTA*xbd%RT2_MS
       UQ_RE = UQ_RE             + 0.5*UQ_RTA*xbd%RT1_RE
     &                           + 0.5*UQ_RTA*xbd%RT2_RE
C
       SCC = bld%SCCON*1.333/(1.0+USA)
       SCC_USA = -SCC/(1.0+USA)
C
       SCC_US1 = SCC_USA*0.5
       SCC_US2 = SCC_USA*0.5
C
C
       SLOG = LOG(xbd%S2/xbd%S1)
       DXI = xbd%X2 - xbd%X1
C
       REZC = SCC*(CQA - SA*ALD)*DXI
     &      - DEA*2.0*          SLOG
     &      + DEA*2.0*(UQ*DXI - ULOG)*bld%DUXCON
C

c        if(  ! (rt2.gt.1.0e3 .and. rt1.le.1.0e3) .or.
c     &     (rt2.gt.1.0e4 .and. rt1.le.1.0e4) .or.
c     &     (rt2.gt.1.0e5 .and. rt1.le.1.0e5)        ) then
c           gga = (HKA-1.0-GCC/RTA)/HKA / sqrt(0.5*CFA)
c           write(*,4455) rta, hka, gga, cfa, cqa, sa, uq, ulog/dxi
c 4455      format(1x,f7.0, 2f9.4,f10.6,2f8.5,2f10.5)
c        endif


       Z_CFA = DEA*2.0*UQ_CFA*DXI * bld%DUXCON
       Z_HKA = DEA*2.0*UQ_HKA*DXI * bld%DUXCON
       Z_DA  = DEA*2.0*UQ_DA *DXI * bld%DUXCON
       Z_SL = -DEA*2.0
       Z_UL = -DEA*2.0 * bld%DUXCON
       Z_DXI = SCC    *(CQA - SA*ALD)     + DEA*2.0*UQ*bld%DUXCON
       Z_USA = SCC_USA*(CQA - SA*ALD)*DXI
       Z_CQA = SCC*DXI
       Z_SA = -SCC*DXI*ALD
       Z_DEA = 2.0*((UQ*DXI - ULOG)*bld%DUXCON - SLOG)
C
       Z_UPW = Z_CQA*(xbd%CQ2-xbd%CQ1) + Z_SA *(xbd%S2 -xbd%S1 )
     &       + Z_CFA*(xbd%CF2-xbd%CF1) + Z_HKA*(xbd%HK2-xbd%HK1)
       Z_DE1 = 0.5*Z_DEA
       Z_DE2 = 0.5*Z_DEA
       Z_US1 = 0.5*Z_USA
       Z_US2 = 0.5*Z_USA
       Z_D1  = 0.5*Z_DA
       Z_D2  = 0.5*Z_DA
       Z_U1  =                 - Z_UL/xbd%U1
       Z_U2  =                   Z_UL/xbd%U2
       Z_X1  = -Z_DXI
       Z_X2  =  Z_DXI
       Z_S1  = (1.0-UPW)*Z_SA  - Z_SL/xbd%S1
       Z_S2  =      UPW *Z_SA  + Z_SL/xbd%S2
       Z_CQ1 = (1.0-UPW)*Z_CQA
       Z_CQ2 =      UPW *Z_CQA
       Z_CF1 = (1.0-UPW)*Z_CFA
       Z_CF2 =      UPW *Z_CFA
       Z_HK1 = (1.0-UPW)*Z_HKA
       Z_HK2 =      UPW *Z_HKA
C
       xbd%VS1(1,1) = Z_S1
       xbd%VS1(1,2) =        Z_UPW*UPW_T1 + Z_DE1*xbd%DE1_T1 + Z_US1
     &  *xbd%US1_T1
       xbd%VS1(1,3) = Z_D1 + Z_UPW*UPW_D1 + Z_DE1*xbd%DE1_D1 + Z_US1
     &  *xbd%US1_D1
       xbd%VS1(1,4) = Z_U1 + Z_UPW*UPW_U1 + Z_DE1*xbd%DE1_U1 + Z_US1
     &  *xbd%US1_U1
       xbd%VS1(1,5) = Z_X1
       xbd%VS2(1,1) = Z_S2
       xbd%VS2(1,2) =        Z_UPW*UPW_T2 + Z_DE2*xbd%DE2_T2 + Z_US2
     &  *xbd%US2_T2
       xbd%VS2(1,3) = Z_D2 + Z_UPW*UPW_D2 + Z_DE2*xbd%DE2_D2 + Z_US2
     &  *xbd%US2_D2
       xbd%VS2(1,4) = Z_U2 + Z_UPW*UPW_U2 + Z_DE2*xbd%DE2_U2 + Z_US2
     &  *xbd%US2_U2
       xbd%VS2(1,5) = Z_X2
       xbd%VSM(1)   =        Z_UPW*UPW_MS + Z_DE1*xbd%DE1_MS + Z_US1
     &  *xbd%US1_MS
     &                                + Z_DE2*xbd%DE2_MS + Z_US2
     &  *xbd%US2_MS
C
       xbd%VS1(1,2) = xbd%VS1(1,2) + Z_CQ1*xbd%CQ1_T1 + Z_CF1*xbd%CF1_T1
     &   + Z_HK1*xbd%HK1_T1
       xbd%VS1(1,3) = xbd%VS1(1,3) + Z_CQ1*xbd%CQ1_D1 + Z_CF1*xbd%CF1_D1
     &   + Z_HK1*xbd%HK1_D1
       xbd%VS1(1,4) = xbd%VS1(1,4) + Z_CQ1*xbd%CQ1_U1 + Z_CF1*xbd%CF1_U1
     &   + Z_HK1*xbd%HK1_U1
C
       xbd%VS2(1,2) = xbd%VS2(1,2) + Z_CQ2*xbd%CQ2_T2 + Z_CF2*xbd%CF2_T2
     &   + Z_HK2*xbd%HK2_T2
       xbd%VS2(1,3) = xbd%VS2(1,3) + Z_CQ2*xbd%CQ2_D2 + Z_CF2*xbd%CF2_D2
     &   + Z_HK2*xbd%HK2_D2
       xbd%VS2(1,4) = xbd%VS2(1,4) + Z_CQ2*xbd%CQ2_U2 + Z_CF2*xbd%CF2_U2
     &   + Z_HK2*xbd%HK2_U2
C
       xbd%VSM(1)   = xbd%VSM(1)   + Z_CQ1*xbd%CQ1_MS + Z_CF1*xbd%CF1_MS
     &   + Z_HK1*xbd%HK1_MS
     &                     + Z_CQ2*xbd%CQ2_MS + Z_CF2*xbd%CF2_MS + Z_HK2
     &  *xbd%HK2_MS
       xbd%VSR(1)   =            Z_CQ1*xbd%CQ1_RE + Z_CF1*xbd%CF1_RE
     &                     + Z_CQ2*xbd%CQ2_RE + Z_CF2*xbd%CF2_RE
       xbd%VSX(1)   = 0.
       xbd%VSREZ(1) = -REZC
C
      ENDIF
C
C**** Set up momentum equation
      HA = 0.5*(xbd%H1 + xbd%H2)
      MA = 0.5*(xbd%M1 + xbd%M2)
      XA = 0.5*(xbd%X1 + xbd%X2)
      TA = 0.5*(xbd%T1 + xbd%T2)
      HWA = 0.5*(xbd%DW1/xbd%T1 + xbd%DW2/xbd%T2)
C
C---- set Cf term, using central value CFM for better accuracy in drag
      CFX     = 0.50*xbd%CFM*XA/TA  +  0.25*(xbd%CF1*xbd%X1/xbd%T1 +
     &   xbd%CF2*xbd%X2/xbd%T2)
      CFX_XA  = 0.50*xbd%CFM   /TA
      CFX_TA  = -.50*xbd%CFM*XA/TA**2
C
      CFX_X1  = 0.25*xbd%CF1   /xbd%T1     + CFX_XA*0.5
      CFX_X2  = 0.25*xbd%CF2   /xbd%T2     + CFX_XA*0.5
      CFX_T1  = -.25*xbd%CF1*xbd%X1/xbd%T1**2  + CFX_TA*0.5
      CFX_T2  = -.25*xbd%CF2*xbd%X2/xbd%T2**2  + CFX_TA*0.5
      CFX_CF1 = 0.25*    xbd%X1/xbd%T1
      CFX_CF2 = 0.25*    xbd%X2/xbd%T2
      CFX_CFM = 0.50*    XA/TA
C
      BTMP = HA + 2.0 - MA + HWA
C
      REZT  = TLOG + BTMP*ULOG - XLOG*0.5*CFX
      Z_CFX = -XLOG*0.5
      Z_HA  =  ULOG
      Z_HWA =  ULOG
      Z_MA  = -ULOG
      Z_XL  =-DDLOG * 0.5*CFX
      Z_UL  = DDLOG * BTMP
      Z_TL  = DDLOG
C
      Z_CFM = Z_CFX*CFX_CFM
      Z_CF1 = Z_CFX*CFX_CF1
      Z_CF2 = Z_CFX*CFX_CF2
C
      Z_T1 = -Z_TL/xbd%T1 + Z_CFX*CFX_T1 + Z_HWA*0.5*(-xbd%DW1/xbd%T1**2
     &  )
      Z_T2 =  Z_TL/xbd%T2 + Z_CFX*CFX_T2 + Z_HWA*0.5*(-xbd%DW2/xbd%T2**2
     &  )
      Z_X1 = -Z_XL/xbd%X1 + Z_CFX*CFX_X1
      Z_X2 =  Z_XL/xbd%X2 + Z_CFX*CFX_X2
      Z_U1 = -Z_UL/xbd%U1
      Z_U2 =  Z_UL/xbd%U2
C
      xbd%VS1(2,2) = 0.5*Z_HA*xbd%H1_T1 + Z_CFM*xbd%CFM_T1 + Z_CF1
     &  *xbd%CF1_T1 + Z_T1
      xbd%VS1(2,3) = 0.5*Z_HA*xbd%H1_D1 + Z_CFM*xbd%CFM_D1 + Z_CF1
     &  *xbd%CF1_D1
      xbd%VS1(2,4) = 0.5*Z_MA*xbd%M1_U1 + Z_CFM*xbd%CFM_U1 + Z_CF1
     &  *xbd%CF1_U1 + Z_U1
      xbd%VS1(2,5) =                                                Z_X1
     &  
      xbd%VS2(2,2) = 0.5*Z_HA*xbd%H2_T2 + Z_CFM*xbd%CFM_T2 + Z_CF2
     &  *xbd%CF2_T2 + Z_T2
      xbd%VS2(2,3) = 0.5*Z_HA*xbd%H2_D2 + Z_CFM*xbd%CFM_D2 + Z_CF2
     &  *xbd%CF2_D2
      xbd%VS2(2,4) = 0.5*Z_MA*xbd%M2_U2 + Z_CFM*xbd%CFM_U2 + Z_CF2
     &  *xbd%CF2_U2 + Z_U2
      xbd%VS2(2,5) =                                                Z_X2
     &  
C
      xbd%VSM(2)   = 0.5*Z_MA*xbd%M1_MS + Z_CFM*xbd%CFM_MS + Z_CF1
     &  *xbd%CF1_MS
     &         + 0.5*Z_MA*xbd%M2_MS                + Z_CF2*xbd%CF2_MS
      xbd%VSR(2)   =                  Z_CFM*xbd%CFM_RE + Z_CF1
     &  *xbd%CF1_RE
     &                                         + Z_CF2*xbd%CF2_RE
      xbd%VSX(2)   = 0.
      xbd%VSREZ(2) = -REZT
C
C**** Set up shape parameter equation
C
      XOT1 = xbd%X1/xbd%T1
      XOT2 = xbd%X2/xbd%T2
C
      HA  = 0.5*(xbd%H1  + xbd%H2 )
      HSA = 0.5*(xbd%HS1 + xbd%HS2)
      HCA = 0.5*(xbd%HC1 + xbd%HC2)
      HWA = 0.5*(xbd%DW1/xbd%T1 + xbd%DW2/xbd%T2)
C
      DIX = (1.0-UPW)*xbd%DI1*XOT1 + UPW*xbd%DI2*XOT2
      CFX = (1.0-UPW)*xbd%CF1*XOT1 + UPW*xbd%CF2*XOT2
      DIX_UPW = xbd%DI2*XOT2 - xbd%DI1*XOT1
      CFX_UPW = xbd%CF2*XOT2 - xbd%CF1*XOT1
C
      BTMP = 2.0*HCA/HSA + 1.0 - HA - HWA
C
      REZH  = HLOG + BTMP*ULOG + XLOG*(0.5*CFX-DIX)
      Z_CFX =  XLOG*0.5
      Z_DIX = -XLOG
      Z_HCA = 2.0*ULOG/HSA
      Z_HA  = -ULOG
      Z_HWA = -ULOG
      Z_XL  = DDLOG * (0.5*CFX-DIX)
      Z_UL  = DDLOG * BTMP
      Z_HL  = DDLOG
C
      Z_UPW = Z_CFX*CFX_UPW + Z_DIX*DIX_UPW
C
      Z_HS1 = -HCA*ULOG/HSA**2 - Z_HL/xbd%HS1
      Z_HS2 = -HCA*ULOG/HSA**2 + Z_HL/xbd%HS2
C
      Z_CF1 = (1.0-UPW)*Z_CFX*XOT1
      Z_CF2 =      UPW *Z_CFX*XOT2
      Z_DI1 = (1.0-UPW)*Z_DIX*XOT1
      Z_DI2 =      UPW *Z_DIX*XOT2
C
      Z_T1 = (1.0-UPW)*(Z_CFX*xbd%CF1 + Z_DIX*xbd%DI1)*(-XOT1/xbd%T1)
      Z_T2 =      UPW *(Z_CFX*xbd%CF2 + Z_DIX*xbd%DI2)*(-XOT2/xbd%T2)
      Z_X1 = (1.0-UPW)*(Z_CFX*xbd%CF1 + Z_DIX*xbd%DI1)/ xbd%T1        -
     &   Z_XL/xbd%X1
      Z_X2 =      UPW *(Z_CFX*xbd%CF2 + Z_DIX*xbd%DI2)/ xbd%T2        +
     &   Z_XL/xbd%X2
      Z_U1 =                                              - Z_UL/xbd%U1
      Z_U2 =                                                Z_UL/xbd%U2
C
      Z_T1 = Z_T1 + Z_HWA*0.5*(-xbd%DW1/xbd%T1**2)
      Z_T2 = Z_T2 + Z_HWA*0.5*(-xbd%DW2/xbd%T2**2)
C
      xbd%VS1(3,1) =                               Z_DI1*xbd%DI1_S1
      xbd%VS1(3,2) = Z_HS1*xbd%HS1_T1 + Z_CF1*xbd%CF1_T1 + Z_DI1
     &  *xbd%DI1_T1 + Z_T1
      xbd%VS1(3,3) = Z_HS1*xbd%HS1_D1 + Z_CF1*xbd%CF1_D1 + Z_DI1
     &  *xbd%DI1_D1
      xbd%VS1(3,4) = Z_HS1*xbd%HS1_U1 + Z_CF1*xbd%CF1_U1 + Z_DI1
     &  *xbd%DI1_U1 + Z_U1
      xbd%VS1(3,5) =                                              Z_X1
      xbd%VS2(3,1) =                               Z_DI2*xbd%DI2_S2
      xbd%VS2(3,2) = Z_HS2*xbd%HS2_T2 + Z_CF2*xbd%CF2_T2 + Z_DI2
     &  *xbd%DI2_T2 + Z_T2
      xbd%VS2(3,3) = Z_HS2*xbd%HS2_D2 + Z_CF2*xbd%CF2_D2 + Z_DI2
     &  *xbd%DI2_D2
      xbd%VS2(3,4) = Z_HS2*xbd%HS2_U2 + Z_CF2*xbd%CF2_U2 + Z_DI2
     &  *xbd%DI2_U2 + Z_U2
      xbd%VS2(3,5) =                                              Z_X2
      xbd%VSM(3)   = Z_HS1*xbd%HS1_MS + Z_CF1*xbd%CF1_MS + Z_DI1
     &  *xbd%DI1_MS
     &         + Z_HS2*xbd%HS2_MS + Z_CF2*xbd%CF2_MS + Z_DI2*xbd%DI2_MS
      xbd%VSR(3)   = Z_HS1*xbd%HS1_RE + Z_CF1*xbd%CF1_RE + Z_DI1
     &  *xbd%DI1_RE
     &         + Z_HS2*xbd%HS2_RE + Z_CF2*xbd%CF2_RE + Z_DI2*xbd%DI2_RE
C
      xbd%VS1(3,2) = xbd%VS1(3,2) + 0.5*(Z_HCA*xbd%HC1_T1+Z_HA*xbd%H1_T1
     &  ) + Z_UPW*UPW_T1
      xbd%VS1(3,3) = xbd%VS1(3,3) + 0.5*(Z_HCA*xbd%HC1_D1+Z_HA*xbd%H1_D1
     &  ) + Z_UPW*UPW_D1
      xbd%VS1(3,4) = xbd%VS1(3,4) + 0.5*(Z_HCA*xbd%HC1_U1           ) +
     &   Z_UPW*UPW_U1
      xbd%VS2(3,2) = xbd%VS2(3,2) + 0.5*(Z_HCA*xbd%HC2_T2+Z_HA*xbd%H2_T2
     &  ) + Z_UPW*UPW_T2
      xbd%VS2(3,3) = xbd%VS2(3,3) + 0.5*(Z_HCA*xbd%HC2_D2+Z_HA*xbd%H2_D2
     &  ) + Z_UPW*UPW_D2
      xbd%VS2(3,4) = xbd%VS2(3,4) + 0.5*(Z_HCA*xbd%HC2_U2           ) +
     &   Z_UPW*UPW_U2
C
      xbd%VSM(3)   = xbd%VSM(3)   + 0.5*(Z_HCA*xbd%HC1_MS           ) +
     &   Z_UPW*UPW_MS
     &                    + 0.5*(Z_HCA*xbd%HC2_MS           )
C
      xbd%VSX(3)   = 0.
      xbd%VSREZ(3) = -REZH
C
      RETURN
      END

C===================================================================70
C
C     Sets up the Newton system governing the
C     transition interval.  Equations governing
C     the  laminar  part  X1 < xi < XT  and
C     the turbulent part  XT < xi < X2
C     are simply summed.
C
C===================================================================70
      SUBROUTINE TRDIF(bld,xbd)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
      REAL*8 :: BL1(4,5), BL2(4,5), BLREZ(4), BLM(4), BLR(4), BLX(4)
     &    , BT1(4,5), BT2(4,5), BTREZ(4), BTM(4), BTR(4), BTX(4)
C
C---- save variables and sensitivities for future restoration
C     DP mod: to remove need for EQUIVALENCE and COM1, COM2
      call store_c1sav(xbd)
      call store_c2sav(xbd)
C      DO 5 ICOM=1, NCOM
C        C1SAV(ICOM) = COM1(ICOM)
C        C2SAV(ICOM) = COM2(ICOM)
C    5 CONTINUE
C
C---- weighting factors for linear interpolation to transition point
      WF2    = (xbd%XT-xbd%X1)/(xbd%X2-xbd%X1)
      WF2_XT = 1.0/(xbd%X2-xbd%X1)
C
      WF2_A1 = WF2_XT*xbd%XT_A1
      WF2_X1 = WF2_XT*xbd%XT_X1 + (WF2-1.0)/(xbd%X2-xbd%X1)
      WF2_X2 = WF2_XT*xbd%XT_X2 -  WF2     /(xbd%X2-xbd%X1)
      WF2_T1 = WF2_XT*xbd%XT_T1
      WF2_T2 = WF2_XT*xbd%XT_T2
      WF2_D1 = WF2_XT*xbd%XT_D1
      WF2_D2 = WF2_XT*xbd%XT_D2
      WF2_U1 = WF2_XT*xbd%XT_U1
      WF2_U2 = WF2_XT*xbd%XT_U2
      WF2_MS = WF2_XT*xbd%XT_MS
      WF2_RE = WF2_XT*xbd%XT_RE
      WF2_XF = WF2_XT*xbd%XT_XF
C
      WF1    = 1.0 - WF2
      WF1_A1 = -WF2_A1
      WF1_X1 = -WF2_X1
      WF1_X2 = -WF2_X2
      WF1_T1 = -WF2_T1
      WF1_T2 = -WF2_T2
      WF1_D1 = -WF2_D1
      WF1_D2 = -WF2_D2
      WF1_U1 = -WF2_U1
      WF1_U2 = -WF2_U2
      WF1_MS = -WF2_MS
      WF1_RE = -WF2_RE
      WF1_XF = -WF2_XF
C
C
C**** FIRST,  do laminar part between X1 and XT
C
C-----interpolate primary variables to transition point
      TT    = xbd%T1*WF1    + xbd%T2*WF2
      TT_A1 = xbd%T1*WF1_A1 + xbd%T2*WF2_A1
      TT_X1 = xbd%T1*WF1_X1 + xbd%T2*WF2_X1
      TT_X2 = xbd%T1*WF1_X2 + xbd%T2*WF2_X2
      TT_T1 = xbd%T1*WF1_T1 + xbd%T2*WF2_T1 + WF1
      TT_T2 = xbd%T1*WF1_T2 + xbd%T2*WF2_T2 + WF2
      TT_D1 = xbd%T1*WF1_D1 + xbd%T2*WF2_D1
      TT_D2 = xbd%T1*WF1_D2 + xbd%T2*WF2_D2
      TT_U1 = xbd%T1*WF1_U1 + xbd%T2*WF2_U1
      TT_U2 = xbd%T1*WF1_U2 + xbd%T2*WF2_U2
      TT_MS = xbd%T1*WF1_MS + xbd%T2*WF2_MS
      TT_RE = xbd%T1*WF1_RE + xbd%T2*WF2_RE
      TT_XF = xbd%T1*WF1_XF + xbd%T2*WF2_XF
C
      DT    = xbd%D1*WF1    + xbd%D2*WF2
      DT_A1 = xbd%D1*WF1_A1 + xbd%D2*WF2_A1
      DT_X1 = xbd%D1*WF1_X1 + xbd%D2*WF2_X1
      DT_X2 = xbd%D1*WF1_X2 + xbd%D2*WF2_X2
      DT_T1 = xbd%D1*WF1_T1 + xbd%D2*WF2_T1
      DT_T2 = xbd%D1*WF1_T2 + xbd%D2*WF2_T2
      DT_D1 = xbd%D1*WF1_D1 + xbd%D2*WF2_D1 + WF1
      DT_D2 = xbd%D1*WF1_D2 + xbd%D2*WF2_D2 + WF2
      DT_U1 = xbd%D1*WF1_U1 + xbd%D2*WF2_U1
      DT_U2 = xbd%D1*WF1_U2 + xbd%D2*WF2_U2
      DT_MS = xbd%D1*WF1_MS + xbd%D2*WF2_MS
      DT_RE = xbd%D1*WF1_RE + xbd%D2*WF2_RE
      DT_XF = xbd%D1*WF1_XF + xbd%D2*WF2_XF
C
      UT    = xbd%U1*WF1    + xbd%U2*WF2
      UT_A1 = xbd%U1*WF1_A1 + xbd%U2*WF2_A1
      UT_X1 = xbd%U1*WF1_X1 + xbd%U2*WF2_X1
      UT_X2 = xbd%U1*WF1_X2 + xbd%U2*WF2_X2
      UT_T1 = xbd%U1*WF1_T1 + xbd%U2*WF2_T1
      UT_T2 = xbd%U1*WF1_T2 + xbd%U2*WF2_T2
      UT_D1 = xbd%U1*WF1_D1 + xbd%U2*WF2_D1
      UT_D2 = xbd%U1*WF1_D2 + xbd%U2*WF2_D2
      UT_U1 = xbd%U1*WF1_U1 + xbd%U2*WF2_U1 + WF1
      UT_U2 = xbd%U1*WF1_U2 + xbd%U2*WF2_U2 + WF2
      UT_MS = xbd%U1*WF1_MS + xbd%U2*WF2_MS
      UT_RE = xbd%U1*WF1_RE + xbd%U2*WF2_RE
      UT_XF = xbd%U1*WF1_XF + xbd%U2*WF2_XF
C
C---- set primary "T" variables at XT  (really placed into "2" variables)
      xbd%X2 = xbd%XT
      xbd%T2 = TT
      xbd%D2 = DT
      xbd%U2 = UT
C
      xbd%AMPL2 = xbd%AMCRIT
      xbd%S2 = 0.
C
C---- calculate laminar secondary "T" variables
      CALL BLKIN(xbd)
      CALL BLVAR(bld,xbd,1)
C
C---- calculate X1-XT midpoint CFM value
      CALL BLMID(bld,xbd,1)
C=
C=    at this point, all "2" variables are really "T" variables at XT
C=
C
C---- set up Newton system for dAm, dTh, dDs, dUe, dXi  at  X1 and XT
      CALL BLDIF(bld,xbd,1)
C
C---- The current Newton system is in terms of "1" and "T" variables,
C-    so calculate its equivalent in terms of "1" and "2" variables.
C-    In other words, convert residual sensitivities wrt "T" variables
C-    into sensitivities wrt "1" and "2" variables.  The amplification
C-    equation is unnecessary here, so the K=1 row is left empty.
      DO 10 K=2, 3
        BLREZ(K) = xbd%VSREZ(K)
        BLM(K)   = xbd%VSM(K)
     &           + xbd%VS2(K,2)*TT_MS
     &           + xbd%VS2(K,3)*DT_MS
     &           + xbd%VS2(K,4)*UT_MS
     &           + xbd%VS2(K,5)*xbd%XT_MS
        BLR(K)   = xbd%VSR(K)
     &           + xbd%VS2(K,2)*TT_RE
     &           + xbd%VS2(K,3)*DT_RE
     &           + xbd%VS2(K,4)*UT_RE
     &           + xbd%VS2(K,5)*xbd%XT_RE
        BLX(K)   = xbd%VSX(K)
     &           + xbd%VS2(K,2)*TT_XF
     &           + xbd%VS2(K,3)*DT_XF
     &           + xbd%VS2(K,4)*UT_XF
     &           + xbd%VS2(K,5)*xbd%XT_XF
C
        BL1(K,1) = xbd%VS1(K,1)
     &           + xbd%VS2(K,2)*TT_A1
     &           + xbd%VS2(K,3)*DT_A1
     &           + xbd%VS2(K,4)*UT_A1
     &           + xbd%VS2(K,5)*xbd%XT_A1
        BL1(K,2) = xbd%VS1(K,2)
     &           + xbd%VS2(K,2)*TT_T1
     &           + xbd%VS2(K,3)*DT_T1
     &           + xbd%VS2(K,4)*UT_T1
     &           + xbd%VS2(K,5)*xbd%XT_T1
        BL1(K,3) = xbd%VS1(K,3)
     &           + xbd%VS2(K,2)*TT_D1
     &           + xbd%VS2(K,3)*DT_D1
     &           + xbd%VS2(K,4)*UT_D1
     &           + xbd%VS2(K,5)*xbd%XT_D1
        BL1(K,4) = xbd%VS1(K,4)
     &           + xbd%VS2(K,2)*TT_U1
     &           + xbd%VS2(K,3)*DT_U1
     &           + xbd%VS2(K,4)*UT_U1
     &           + xbd%VS2(K,5)*xbd%XT_U1
        BL1(K,5) = xbd%VS1(K,5)
     &           + xbd%VS2(K,2)*TT_X1
     &           + xbd%VS2(K,3)*DT_X1
     &           + xbd%VS2(K,4)*UT_X1
     &           + xbd%VS2(K,5)*xbd%XT_X1
C
        BL2(K,1) = 0.
        BL2(K,2) = xbd%VS2(K,2)*TT_T2
     &           + xbd%VS2(K,3)*DT_T2
     &           + xbd%VS2(K,4)*UT_T2
     &           + xbd%VS2(K,5)*xbd%XT_T2
        BL2(K,3) = xbd%VS2(K,2)*TT_D2
     &           + xbd%VS2(K,3)*DT_D2
     &           + xbd%VS2(K,4)*UT_D2
     &           + xbd%VS2(K,5)*xbd%XT_D2
        BL2(K,4) = xbd%VS2(K,2)*TT_U2
     &           + xbd%VS2(K,3)*DT_U2
     &           + xbd%VS2(K,4)*UT_U2
     &           + xbd%VS2(K,5)*xbd%XT_U2
        BL2(K,5) = xbd%VS2(K,2)*TT_X2
     &           + xbd%VS2(K,3)*DT_X2
     &           + xbd%VS2(K,4)*UT_X2
     &           + xbd%VS2(K,5)*xbd%XT_X2
C
   10 CONTINUE
C
C
C**** SECOND, set up turbulent part between XT and X2  ****
C
C---- calculate equilibrium shear coefficient CQT at transition point
      CALL BLVAR(bld,xbd,2)
C
C---- set initial shear coefficient value ST at transition point
C-    ( note that CQ2, CQ2_T2, etc. are really "CQT", "CQT_TT", etc.)
C
      CTR     = bld%CTRCON*EXP(-bld%CTRCEX/(xbd%HK2-1.0))
      CTR_HK2 = CTR * bld%CTRCEX/(xbd%HK2-1.0)**2
C
c      CTR     = 1.1*EXP(-10.0/HK2**2)
c      CTR_HK2 = CTR * 10.0 * 2.0/HK2**3
C
CCC      CTR = 1.2
CCC      CTR = 0.7
CCC      CTR_HK2 = 0.0
C
      ST    = CTR*xbd%CQ2
      ST_TT = CTR*xbd%CQ2_T2 + xbd%CQ2*CTR_HK2*xbd%HK2_T2
      ST_DT = CTR*xbd%CQ2_D2 + xbd%CQ2*CTR_HK2*xbd%HK2_D2
      ST_UT = CTR*xbd%CQ2_U2 + xbd%CQ2*CTR_HK2*xbd%HK2_U2
      ST_MS = CTR*xbd%CQ2_MS + xbd%CQ2*CTR_HK2*xbd%HK2_MS
      ST_RE = CTR*xbd%CQ2_RE
C
C---- calculate ST sensitivities wrt the actual "1" and "2" variables
      ST_A1 = ST_TT*TT_A1 + ST_DT*DT_A1 + ST_UT*UT_A1
      ST_X1 = ST_TT*TT_X1 + ST_DT*DT_X1 + ST_UT*UT_X1
      ST_X2 = ST_TT*TT_X2 + ST_DT*DT_X2 + ST_UT*UT_X2
      ST_T1 = ST_TT*TT_T1 + ST_DT*DT_T1 + ST_UT*UT_T1
      ST_T2 = ST_TT*TT_T2 + ST_DT*DT_T2 + ST_UT*UT_T2
      ST_D1 = ST_TT*TT_D1 + ST_DT*DT_D1 + ST_UT*UT_D1
      ST_D2 = ST_TT*TT_D2 + ST_DT*DT_D2 + ST_UT*UT_D2
      ST_U1 = ST_TT*TT_U1 + ST_DT*DT_U1 + ST_UT*UT_U1
      ST_U2 = ST_TT*TT_U2 + ST_DT*DT_U2 + ST_UT*UT_U2
      ST_MS = ST_TT*TT_MS + ST_DT*DT_MS + ST_UT*UT_MS + ST_MS
      ST_RE = ST_TT*TT_RE + ST_DT*DT_RE + ST_UT*UT_RE + ST_RE
      ST_XF = ST_TT*TT_XF + ST_DT*DT_XF + ST_UT*UT_XF
C
      xbd%AMPL2 = 0.
      xbd%S2 = ST
C
C---- recalculate turbulent secondary "T" variables using proper CTI
      CALL BLVAR(bld,xbd,2)
C
C---- set "1" variables to "T" variables and reset "2" variables
C-    to their saved turbulent values
C
C     DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
      call com2_to_com1(xbd)
      call from_c2sav(xbd)
C      DO 30 ICOM=1, NCOM
C        COM1(ICOM) = COM2(ICOM)
C        COM2(ICOM) = C2SAV(ICOM)
C   30 CONTINUE
C
C---- calculate XT-X2 midpoint CFM value
      CALL BLMID(bld,xbd,2)
C
C---- set up Newton system for dCt, dTh, dDs, dUe, dXi  at  XT and X2
      CALL BLDIF(bld,xbd,2)
C
C---- convert sensitivities wrt "T" variables into sensitivities
C-    wrt "1" and "2" variables as done before for the laminar part
      DO 40 K=1, 3
        BTREZ(K) = xbd%VSREZ(K)
        BTM(K)   = xbd%VSM(K) 
     &           + xbd%VS1(K,1)*ST_MS
     &           + xbd%VS1(K,2)*TT_MS
     &           + xbd%VS1(K,3)*DT_MS
     &           + xbd%VS1(K,4)*UT_MS
     &           + xbd%VS1(K,5)*xbd%XT_MS
        BTR(K)   = xbd%VSR(K) 
     &           + xbd%VS1(K,1)*ST_RE
     &           + xbd%VS1(K,2)*TT_RE
     &           + xbd%VS1(K,3)*DT_RE
     &           + xbd%VS1(K,4)*UT_RE
     &           + xbd%VS1(K,5)*xbd%XT_RE
        BTX(K)   = xbd%VSX(K)
     &           + xbd%VS1(K,1)*ST_XF
     &           + xbd%VS1(K,2)*TT_XF
     &           + xbd%VS1(K,3)*DT_XF
     &           + xbd%VS1(K,4)*UT_XF
     &           + xbd%VS1(K,5)*xbd%XT_XF
C
        BT1(K,1) = xbd%VS1(K,1)*ST_A1
     &           + xbd%VS1(K,2)*TT_A1
     &           + xbd%VS1(K,3)*DT_A1
     &           + xbd%VS1(K,4)*UT_A1
     &           + xbd%VS1(K,5)*xbd%XT_A1
        BT1(K,2) = xbd%VS1(K,1)*ST_T1
     &           + xbd%VS1(K,2)*TT_T1
     &           + xbd%VS1(K,3)*DT_T1
     &           + xbd%VS1(K,4)*UT_T1
     &           + xbd%VS1(K,5)*xbd%XT_T1
        BT1(K,3) = xbd%VS1(K,1)*ST_D1
     &           + xbd%VS1(K,2)*TT_D1
     &           + xbd%VS1(K,3)*DT_D1
     &           + xbd%VS1(K,4)*UT_D1
     &           + xbd%VS1(K,5)*xbd%XT_D1
        BT1(K,4) = xbd%VS1(K,1)*ST_U1
     &           + xbd%VS1(K,2)*TT_U1
     &           + xbd%VS1(K,3)*DT_U1
     &           + xbd%VS1(K,4)*UT_U1
     &           + xbd%VS1(K,5)*xbd%XT_U1
        BT1(K,5) = xbd%VS1(K,1)*ST_X1
     &           + xbd%VS1(K,2)*TT_X1
     &           + xbd%VS1(K,3)*DT_X1
     &           + xbd%VS1(K,4)*UT_X1
     &           + xbd%VS1(K,5)*xbd%XT_X1
C
        BT2(K,1) = xbd%VS2(K,1)
        BT2(K,2) = xbd%VS2(K,2)
     &           + xbd%VS1(K,1)*ST_T2
     &           + xbd%VS1(K,2)*TT_T2
     &           + xbd%VS1(K,3)*DT_T2
     &           + xbd%VS1(K,4)*UT_T2
     &           + xbd%VS1(K,5)*xbd%XT_T2
        BT2(K,3) = xbd%VS2(K,3)
     &           + xbd%VS1(K,1)*ST_D2
     &           + xbd%VS1(K,2)*TT_D2
     &           + xbd%VS1(K,3)*DT_D2
     &           + xbd%VS1(K,4)*UT_D2
     &           + xbd%VS1(K,5)*xbd%XT_D2
        BT2(K,4) = xbd%VS2(K,4)
     &           + xbd%VS1(K,1)*ST_U2
     &           + xbd%VS1(K,2)*TT_U2
     &           + xbd%VS1(K,3)*DT_U2
     &           + xbd%VS1(K,4)*UT_U2
     &           + xbd%VS1(K,5)*xbd%XT_U2
        BT2(K,5) = xbd%VS2(K,5)
     &           + xbd%VS1(K,1)*ST_X2
     &           + xbd%VS1(K,2)*TT_X2
     &           + xbd%VS1(K,3)*DT_X2
     &           + xbd%VS1(K,4)*UT_X2
     &           + xbd%VS1(K,5)*xbd%XT_X2
C
   40 CONTINUE
C
C---- Add up laminar and turbulent parts to get final system
C-    in terms of honest-to-God "1" and "2" variables.
      xbd%VSREZ(1) =            BTREZ(1)
      xbd%VSREZ(2) = BLREZ(2) + BTREZ(2)
      xbd%VSREZ(3) = BLREZ(3) + BTREZ(3)
      xbd%VSM(1)   =            BTM(1)
      xbd%VSM(2)   = BLM(2)   + BTM(2)
      xbd%VSM(3)   = BLM(3)   + BTM(3)
      xbd%VSR(1)   =            BTR(1)
      xbd%VSR(2)   = BLR(2)   + BTR(2)
      xbd%VSR(3)   = BLR(3)   + BTR(3)
      xbd%VSX(1)   =            BTX(1)
      xbd%VSX(2)   = BLX(2)   + BTX(2)
      xbd%VSX(3)   = BLX(3)   + BTX(3)
      DO 60 L=1, 5
        xbd%VS1(1,L) =            BT1(1,L)
        xbd%VS2(1,L) =            BT2(1,L)
        xbd%VS1(2,L) = BL1(2,L) + BT1(2,L)
        xbd%VS2(2,L) = BL2(2,L) + BT2(2,L)
        xbd%VS1(3,L) = BL1(3,L) + BT1(3,L)
        xbd%VS2(3,L) = BL2(3,L) + BT2(3,L)
   60 CONTINUE
C
C---- To be sanitary, restore "1" quantities which got clobbered
C-    in all of the numerical gymnastics above.  The "2" variables
C-    were already restored for the XT-X2 differencing part.
C     DP mod: to remove need for EQUIVALENCE and COM1
      call from_c1sav(xbd)
C      DO 70 ICOM=1, NCOM
C        COM1(ICOM) = C1SAV(ICOM)
C   70 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C
C     Sets up the BL Newton system governing the current interval:
C
C     |       ||dA1|     |       ||dA2|       |     |
C     |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
C     |       ||dD1|     |       ||dD2|       |     |
C              |dU1|              |dU2|
C              |dX1|              |dX2|
C
C        3x5    5x1         3x5    5x1          3x1
C
C     The system as shown corresponds to a laminar station
C     If TRAN, then  dS2  replaces  dA2
C     If TURB, then  dS1, dS2  replace  dA1, dA2
C
C
C===================================================================70
      SUBROUTINE BLSYS(bld,xbd)

      use xfoil_data_mod
      IMPLICIT REAL(M)
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
C---- calculate secondary BL variables and their sensitivities
      IF(xbd%WAKE) THEN
       CALL BLVAR(bld,xbd,3)
       CALL BLMID(bld,xbd,3)
      ELSE IF(xbd%TURB.OR.xbd%TRAN) THEN
       CALL BLVAR(bld,xbd,2)
       CALL BLMID(bld,xbd,2)
      ELSE
       CALL BLVAR(bld,xbd,1)
       CALL BLMID(bld,xbd,1)
      ENDIF
C
C---- for the similarity station, "1" and "2" variables are the same
      IF(xbd%SIMI) THEN
C      DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
       call com2_to_com1(xbd)
C       DO 3 ICOM=1, NCOM
C         COM1(ICOM) = COM2(ICOM)
C    3  CONTINUE
      ENDIF
C
C---- set up appropriate finite difference system for current interval
      IF(xbd%TRAN) THEN
       CALL TRDIF(bld,xbd)
      ELSE IF(xbd%SIMI) THEN
       CALL BLDIF(bld,xbd,0)
      ELSE IF(.NOT.xbd%TURB) THEN
       CALL BLDIF(bld,xbd,1)
      ELSE IF(xbd%WAKE) THEN
       CALL BLDIF(bld,xbd,3)
      ELSE IF(xbd%TURB) THEN
       CALL BLDIF(bld,xbd,2)
      ENDIF
C
      IF(xbd%SIMI) THEN
C----- at similarity station, "1" variables are really "2" variables
       DO 10 K=1, 4
         DO 101 L=1, 5
           xbd%VS2(K,L) = xbd%VS1(K,L) + xbd%VS2(K,L)
           xbd%VS1(K,L) = 0.
  101    CONTINUE
   10  CONTINUE
      ENDIF
C
C---- change system over into incompressible Uei and Mach
      DO 20 K=1, 4
C
C------ residual derivatives wrt compressible Uec
        RES_U1 = xbd%VS1(K,4)
        RES_U2 = xbd%VS2(K,4)
        RES_MS = xbd%VSM(K)
C
C------ combine with derivatives of compressible  U1,U2 = Uec(Uei M)
        xbd%VS1(K,4) = RES_U1*xbd%U1_UEI
        xbd%VS2(K,4) =                RES_U2*xbd%U2_UEI
        xbd%VSM(K)   = RES_U1*xbd%U1_MS + RES_U2*xbd%U2_MS  + RES_MS
   20 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE DSLIM(DSTR,THET,UEDG,MSQ,HKLIM)
      IMPLICIT REAL (A-H,M,O-Z)
C
      H = DSTR/THET
      CALL HKIN(H,MSQ,HK,HK_H,HK_M)
C
      DH = MAX( 0.0 , HKLIM-HK ) / HK_H
      DSTR = DSTR + DH*THET
C
      RETURN
      END

C===================================================================70
C
C     Marches the BLs and wake in direct mode using
C     the UEDG array. If separation is encountered,
C     a plausible value of Hk extrapolated from
C     upstream is prescribed instead.  Continuous
C     checking of transition onset is performed.
C
C===================================================================70
      SUBROUTINE MRCHUE(xfd,bld,xbd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd

      LOGICAL DIRECT
      REAL*8 MSQ
C
C---- shape parameters for separation criteria
      HLMAX = 3.8
      HTMAX = 2.5
C
      DO 2000 IS=1, 2
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE) WRITE(*,*) '   side ', IS, ' ...'

      xbd%AMCRIT = xfd%ACRIT
C
C---- set forced transition arc length position
      CALL XIFSET(xfd,xbd,IS)
C
C---- initialize similarity station with Thwaites' formula
      IBL = 2
      XSI = xfd%XSSI(IBL,IS)
      UEI = xfd%UEDG(IBL,IS)
C      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
C      BULE = MAX( -.08 , BULE )
      xbd%BULE = 1.0
      UCON = UEI/XSI**xbd%BULE
      TSQ = 0.45/(UCON*(5.0*xbd%BULE+1.0)*xbd%REYBL) * XSI**(1.0
     &  -xbd%BULE)
      THI = SQRT(TSQ)
      DSI = 2.2*THI
      AMI = 0.0
C
C---- initialize Ctau for first turbulent station
      CTI = 0.03
C
      xbd%TRAN = .FALSE.
      xbd%TURB = .FALSE.
      xfd%ITRAN(IS) = xfd%IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, xfd%NBL(IS)
        IBM = IBL-1
C
        IW = IBL - xfd%IBLTE(IS)
C
        xbd%SIMI = IBL.EQ.2
        xbd%WAKE = IBL.GT.xfd%IBLTE(IS)
C
C------ prescribed quantities
        XSI = xfd%XSSI(IBL,IS)
        UEI = xfd%UEDG(IBL,IS)
C
        IF(xbd%WAKE) THEN
         IW = IBL - xfd%IBLTE(IS)
         DSWAKI = xfd%WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        DIRECT = .TRUE.
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
C
          CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN(xbd)
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.xbd%SIMI) .AND. (.NOT.xbd%TURB)) THEN
           CALL TRCHEK(xfd,xbd)
C          DP mod: check for infinite loop condition
           IF (xfd%XFOIL_FAIL) RETURN
           AMI = xbd%AMPL2
C
           IF(xbd%TRAN) THEN
            xfd%ITRAN(IS) = IBL
            IF(CTI.LE.0.0) THEN
             CTI = 0.03
             xbd%S2 = CTI
            ENDIF
           ELSE
            xfd%ITRAN(IS) = IBL+2
           ENDIF
C
C
          ENDIF
C
          IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
           TTE = xfd%THET(xfd%IBLTE(1),1) + xfd%THET(xfd%IBLTE(2),2)
           DTE = xfd%DSTR(xfd%IBLTE(1),1) + xfd%DSTR(xfd%IBLTE(2),2) +
     &   xfd%ANTE
           CTE = ( xfd%CTAU(xfd%IBLTE(1),1)*xfd%THET(xfd%IBLTE(1),1)
     &           + xfd%CTAU(xfd%IBLTE(2),2)*xfd%THET(xfd%IBLTE(2),2) ) /
     &   TTE
           CALL TESYS(bld,xbd,CTE,TTE,DTE)
          ELSE
           CALL BLSYS(bld,xbd)
          ENDIF
C
          IF(DIRECT) THEN
C
C--------- try direct mode (set dUe = 0 in currently empty 4th line)
           xbd%VS2(4,1) = 0.
           xbd%VS2(4,2) = 0.
           xbd%VS2(4,3) = 0.
           xbd%VS2(4,4) = 1.0
           xbd%VSREZ(4) = 0.
C
C--------- solve Newton system for current "2" station
           CALL GAUSS(4,4,xbd%VS2,xbd%VSREZ,1)
C
C--------- determine max changes and underrelax if necessary
           DMAX = MAX( ABS(xbd%VSREZ(2)/THI),
     &                 ABS(xbd%VSREZ(3)/DSI)  )
           IF(IBL.LT.xfd%ITRAN(IS)) DMAX = MAX(DMAX,ABS(xbd%VSREZ(1)/10
     &  .0))
           IF(IBL.GE.xfd%ITRAN(IS)) DMAX = MAX(DMAX,ABS(xbd%VSREZ(1)/CTI
     &   ))
C
           xfd%RLX = 1.0
           IF(DMAX.GT.0.3) xfd%RLX = 0.3/DMAX
C
C--------- see if direct mode is not applicable
           IF(IBL .NE. xfd%IBLTE(IS)+1) THEN
C
C---------- calculate resulting kinematic shape parameter Hk
            MSQ = UEI*UEI*xbd%HSTINV / (xbd%GM1BL*(1.0 - 0.5*UEI*UEI
     &  *xbd%HSTINV))
            HTEST = (DSI + xfd%RLX*xbd%VSREZ(3)) / (THI + xfd%RLX
     &  *xbd%VSREZ(2))
            CALL HKIN( HTEST, MSQ, HKTEST, DUMMY, DUMMY)
C
C---------- decide whether to do direct or inverse problem based on Hk
            IF(IBL.LT.xfd%ITRAN(IS)) HMAX = HLMAX
            IF(IBL.GE.xfd%ITRAN(IS)) HMAX = HTMAX
            DIRECT = HKTEST.LT.HMAX
           ENDIF
C
           IF(DIRECT) THEN
C---------- update as usual
ccc            IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
            IF(IBL.GE.xfd%ITRAN(IS)) CTI = CTI + xfd%RLX*xbd%VSREZ(1)
            THI = THI + xfd%RLX*xbd%VSREZ(2)
            DSI = DSI + xfd%RLX*xbd%VSREZ(3)
           ELSE
C---------- set prescribed Hk for inverse calculation at the current station
            IF(IBL.LT.xfd%ITRAN(IS)) THEN
C----------- laminar case: relatively slow increase in Hk downstream
             HTARG = xbd%HK1 + 0.03*(xbd%X2-xbd%X1)/xbd%T1
            ELSE IF(IBL.EQ.xfd%ITRAN(IS)) THEN
C----------- transition interval: weighted laminar and turbulent case
             HTARG = xbd%HK1 + (0.03*(xbd%XT-xbd%X1) - 0.15*(xbd%X2
     &  -xbd%XT))/xbd%T1
            ELSE IF(xbd%WAKE) THEN
C----------- turbulent wake case:
C-           asymptotic wake behavior with approximate Backward Euler
             CONST = 0.03*(xbd%X2-xbd%X1)/xbd%T1
             xbd%HK2 = xbd%HK1
             xbd%HK2 = xbd%HK2 - (xbd%HK2 +     CONST*(xbd%HK2-1.0)**3 -
     &   xbd%HK1)
     &                  /(1.0 + 3.0*CONST*(xbd%HK2-1.0)**2)
             xbd%HK2 = xbd%HK2 - (xbd%HK2 +     CONST*(xbd%HK2-1.0)**3 -
     &   xbd%HK1)
     &                  /(1.0 + 3.0*CONST*(xbd%HK2-1.0)**2)
             xbd%HK2 = xbd%HK2 - (xbd%HK2 +     CONST*(xbd%HK2-1.0)**3 -
     &   xbd%HK1)
     &                  /(1.0 + 3.0*CONST*(xbd%HK2-1.0)**2)
             HTARG = xbd%HK2
            ELSE
C----------- turbulent case: relatively fast decrease in Hk downstream
             HTARG = xbd%HK1 - 0.15*(xbd%X2-xbd%X1)/xbd%T1
            ENDIF
C
C---------- limit specified Hk to something reasonable
            IF(xbd%WAKE) THEN
             HTARG = MAX( HTARG , 1.01 )
            ELSE
             HTARG = MAX( HTARG , HMAX )
            ENDIF
C
C           DP mod: added SILENT_MODE option
C           DP mod: change write precision
            IF (.NOT. xfd%SILENT_MODE) WRITE(*,1300) IBL, HTARG
 1300       FORMAT(' MRCHUE: Inverse mode at', I4, '     Hk =', F15.12)
C
C---------- try again with prescribed Hk
            GO TO 100
C
           ENDIF
C
          ELSE
C
C-------- inverse mode (force Hk to prescribed value HTARG)
           xbd%VS2(4,1) = 0.
           xbd%VS2(4,2) = xbd%HK2_T2
           xbd%VS2(4,3) = xbd%HK2_D2
           xbd%VS2(4,4) = xbd%HK2_U2
           xbd%VSREZ(4) = HTARG - xbd%HK2
C
           CALL GAUSS(4,4,xbd%VS2,xbd%VSREZ,1)
C
C--------- added Ue clamp   MD  3 Apr 03
           DMAX = MAX( ABS(xbd%VSREZ(2)/THI),
     &                 ABS(xbd%VSREZ(3)/DSI),
     &                 ABS(xbd%VSREZ(4)/UEI)  )
           IF(IBL.GE.xfd%ITRAN(IS)) DMAX = MAX( DMAX , ABS(xbd%VSREZ(1)
     &  /CTI))
C
           xfd%RLX = 1.0
           IF(DMAX.GT.0.3) xfd%RLX = 0.3/DMAX
C
C--------- update variables
ccc           IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
           IF(IBL.GE.xfd%ITRAN(IS)) CTI = CTI + xfd%RLX*xbd%VSREZ(1)
           THI = THI + xfd%RLX*xbd%VSREZ(2)
           DSI = DSI + xfd%RLX*xbd%VSREZ(3)
           UEI = UEI + xfd%RLX*xbd%VSREZ(4)
C
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.xfd%ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.xfd%IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*xbd%HSTINV / (xbd%GM1BL*(1.0 - 0.5*UEI*UEI
     &  *xbd%HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.1.0E-5) GO TO 110
C
  100   CONTINUE
C       DP mod: added SILENT_MODE option
        IF (.NOT. xfd%SILENT_MODE) WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHUE: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN 
          IF(IBL.LE.xfd%IBLTE(IS)) THEN
           THI = xfd%THET(IBM,IS) * (xfd%XSSI(IBL,IS)/xfd%XSSI(IBM,IS))*
     &  *0.5
           DSI = xfd%DSTR(IBM,IS) * (xfd%XSSI(IBL,IS)/xfd%XSSI(IBM,IS))*
     &  *0.5
          ELSE IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
          ELSE
           THI = xfd%THET(IBM,IS)
           RATLEN = (xfd%XSSI(IBL,IS)-xfd%XSSI(IBM,IS)) / (10.0
     &  *xfd%DSTR(IBM,IS))
           DSI = (xfd%DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
          ENDIF
          IF(IBL.EQ.xfd%ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.xfd%ITRAN(IS)) CTI = xfd%CTAU(IBM,IS)
C
          UEI = xfd%UEDG(IBL,IS)
          IF(IBL.GT.2 .AND. IBL.LT.xfd%NBL(IS))
     &     UEI = 0.5*(xfd%UEDG(IBL-1,IS) + xfd%UEDG(IBL+1,IS))
         ENDIF
C
 109     CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN(xbd)
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.xbd%SIMI) .AND. (.NOT.xbd%TURB)) THEN
          CALL TRCHEK(xfd,xbd)
C         DP mod: check for infinite loop condition
          IF (xfd%XFOIL_FAIL) RETURN
          AMI = xbd%AMPL2
          IF(     xbd%TRAN) xfd%ITRAN(IS) = IBL
          IF(.NOT.xbd%TRAN) xfd%ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.xfd%ITRAN(IS)) CALL BLVAR(bld,xbd,1)
         IF(IBL.GE.xfd%ITRAN(IS)) CALL BLVAR(bld,xbd,2)
         IF(xbd%WAKE) CALL BLVAR(bld,xbd,3)
C
         IF(IBL.LT.xfd%ITRAN(IS)) CALL BLMID(bld,xbd,1)
         IF(IBL.GE.xfd%ITRAN(IS)) CALL BLMID(bld,xbd,2)
         IF(xbd%WAKE) CALL BLMID(bld,xbd,3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
C------ store primary variables
        IF(IBL.LT.xfd%ITRAN(IS)) xfd%CTAU(IBL,IS) = AMI
        IF(IBL.GE.xfd%ITRAN(IS)) xfd%CTAU(IBL,IS) = CTI
        xfd%THET(IBL,IS) = THI
        xfd%DSTR(IBL,IS) = DSI
        xfd%UEDG(IBL,IS) = UEI
        xfd%MASS(IBL,IS) = DSI*UEI
        xfd%TAU(IBL,IS)  = 0.5*xbd%R2*xbd%U2*xbd%U2*xbd%CF2
        xfd%DIS(IBL,IS)  =     xbd%R2*xbd%U2*xbd%U2*xbd%U2*xbd%DI2
     &  *xbd%HS2*0.5
        xfd%CTQ(IBL,IS)  = xbd%CQ2
        xfd%DELT(IBL,IS) = xbd%DE2
        xfd%TSTR(IBL,IS) = xbd%HS2*xbd%T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN(xbd)
C       DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
        call com2_to_com1(xbd)
C        DO 310 ICOM=1, NCOM
C          COM1(ICOM) = COM2(ICOM)
C  310   CONTINUE
C
C------ turbulent intervals will follow transition interval or TE
        IF(xbd%TRAN .OR. IBL.EQ.xfd%IBLTE(IS)) THEN
         xbd%TURB = .TRUE.
C
C------- save transition location
         xfd%TFORCE(IS) = xbd%TRFORC
         xfd%XSSITR(IS) = xbd%XT
        ENDIF
C
        xbd%TRAN = .FALSE.
C
        IF(IBL.EQ.xfd%IBLTE(IS)) THEN
         THI = xfd%THET(xfd%IBLTE(1),1) + xfd%THET(xfd%IBLTE(2),2)
         DSI = xfd%DSTR(xfd%IBLTE(1),1) + xfd%DSTR(xfd%IBLTE(2),2) +
     &   xfd%ANTE
        ENDIF
C
 1000 CONTINUE
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Marches the BLs and wake in mixed mode using
C     the current Ue and Hk.  The calculated Ue
C     and Hk lie along a line quasi-normal to the
C     natural Ue-Hk characteristic line of the
C     current BL so that the Goldstein or Levy-Lees
C     singularity is never encountered.  Continuous
C     checking of transition onset is performed.
C
C===================================================================70
      SUBROUTINE MRCHDU(xfd,bld,xbd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd

      REAL*8 :: VTMP(4,5), VZTMP(4)
      REAL*8 MSQ
ccc   REAL MDI
C
      DATA DEPS / 5.0E-6 /
C
C---- constant controlling how far Hk is allowed to deviate
C-    from the specified value.
      SENSWT = 1000.0
C
      DO 2000 IS=1, 2

      xbd%AMCRIT = xfd%ACRIT
C
C---- set forced transition arc length position
      CALL XIFSET(xfd,xbd,IS)
C
C---- set leading edge pressure gradient parameter  x/u du/dx
      IBL = 2
      XSI = xfd%XSSI(IBL,IS)
      UEI = xfd%UEDG(IBL,IS)
CCC      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
CCC      BULE = MAX( -.08 , BULE )
      xbd%BULE = 1.0
C
C---- old transition station
      ITROLD = xfd%ITRAN(IS)
C
      xbd%TRAN = .FALSE.
      xbd%TURB = .FALSE.
      xfd%ITRAN(IS) = xfd%IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, xfd%NBL(IS)
        IBM = IBL-1
C
        xbd%SIMI = IBL.EQ.2
        xbd%WAKE = IBL.GT.xfd%IBLTE(IS)
C
C------ initialize current station to existing variables
        XSI = xfd%XSSI(IBL,IS)
        UEI = xfd%UEDG(IBL,IS)
        THI = xfd%THET(IBL,IS)
        DSI = xfd%DSTR(IBL,IS)

CCC        MDI = MASS(IBL,IS)
C
C------ fixed BUG   MD 7 June 99
        IF(IBL.LT.ITROLD) THEN
         AMI = xfd%CTAU(IBL,IS)
         CTI = 0.03
        ELSE
         CTI = xfd%CTAU(IBL,IS)
         IF(CTI.LE.0.0) CTI = 0.03
        ENDIF
C
CCC        DSI = MDI/UEI
C
        IF(xbd%WAKE) THEN
         IW = IBL - xfd%IBLTE(IS)
         DSWAKI = xfd%WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        IF(IBL.LE.xfd%IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.02000*THI) +
     &   DSWAKI
        IF(IBL.GT.xfd%IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.00005*THI) +
     &   DSWAKI
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
          CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN(xbd)
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.xbd%SIMI) .AND. (.NOT.xbd%TURB)) THEN
           CALL TRCHEK(xfd,xbd)
C          DP mod: check for infinite loop condition
           IF (xfd%XFOIL_FAIL) RETURN
           AMI = xbd%AMPL2
           IF(     xbd%TRAN) xfd%ITRAN(IS) = IBL
           IF(.NOT.xbd%TRAN) xfd%ITRAN(IS) = IBL+2
          ENDIF
C
          IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
           TTE = xfd%THET(xfd%IBLTE(1),1) + xfd%THET(xfd%IBLTE(2),2)
           DTE = xfd%DSTR(xfd%IBLTE(1),1) + xfd%DSTR(xfd%IBLTE(2),2) +
     &   xfd%ANTE
           CTE = ( xfd%CTAU(xfd%IBLTE(1),1)*xfd%THET(xfd%IBLTE(1),1)
     &           + xfd%CTAU(xfd%IBLTE(2),2)*xfd%THET(xfd%IBLTE(2),2) ) /
     &   TTE
           CALL TESYS(bld,xbd,CTE,TTE,DTE)
          ELSE
           CALL BLSYS(bld,xbd)
          ENDIF
C
C-------- set stuff at first iteration...
          IF(ITBL.EQ.1) THEN
C
C--------- set "baseline" Ue and Hk for forming  Ue(Hk)  relation
           UEREF = xbd%U2
           HKREF = xbd%HK2
C
C--------- if current point IBL was turbulent and is now laminar, then...
           IF(IBL.LT.xfd%ITRAN(IS) .AND. IBL.GE.ITROLD ) THEN
C---------- extrapolate baseline Hk
            UEM = xfd%UEDG(IBL-1,IS)
            DSM = xfd%DSTR(IBL-1,IS)
            THM = xfd%THET(IBL-1,IS)
            MSQ = UEM*UEM*xbd%HSTINV / (xbd%GM1BL*(1.0 - 0.5*UEM*UEM
     &  *xbd%HSTINV))
            CALL HKIN( DSM/THM, MSQ, HKREF, DUMMY, DUMMY )
           ENDIF
C
C--------- if current point IBL was laminar, then...
           IF(IBL.LT.ITROLD) THEN
C---------- reinitialize or extrapolate Ctau if it's now turbulent
            IF(xbd%TRAN) xfd%CTAU(IBL,IS) = 0.03
            IF(xbd%TURB) xfd%CTAU(IBL,IS) = xfd%CTAU(IBL-1,IS)
            IF(xbd%TRAN .OR. xbd%TURB) THEN
             CTI = xfd%CTAU(IBL,IS)
             xbd%S2 = CTI
            ENDIF
           ENDIF
C
          ENDIF
C
C
          IF(xbd%SIMI .OR. IBL.EQ.xfd%IBLTE(IS)+1) THEN
C
C--------- for similarity station or first wake point, prescribe Ue
           xbd%VS2(4,1) = 0.
           xbd%VS2(4,2) = 0.
           xbd%VS2(4,3) = 0.
           xbd%VS2(4,4) = xbd%U2_UEI
           xbd%VSREZ(4) = UEREF - xbd%U2
C
          ELSE
C
C********* calculate Ue-Hk characteristic slope
C
           DO 20 K=1, 4
             VZTMP(K) = xbd%VSREZ(K)
             DO 201 L=1, 5
               VTMP(K,L) = xbd%VS2(K,L)
  201        CONTINUE
   20      CONTINUE
C
C--------- set unit dHk
           VTMP(4,1) = 0.
           VTMP(4,2) = xbd%HK2_T2
           VTMP(4,3) = xbd%HK2_D2
           VTMP(4,4) = xbd%HK2_U2*xbd%U2_UEI
           VZTMP(4)  = 1.0
C
C--------- calculate dUe response
           CALL GAUSS(4,4,VTMP,VZTMP,1)
C
C--------- set  SENSWT * (normalized dUe/dHk)
           SENNEW = SENSWT * VZTMP(4) * HKREF/UEREF
           IF(ITBL.LE.5) THEN
            SENS = SENNEW
           ELSE IF(ITBL.LE.15) THEN
            SENS = 0.5*(SENS + SENNEW)
           ENDIF
C
C--------- set prescribed Ue-Hk combination
           xbd%VS2(4,1) = 0.
           xbd%VS2(4,2) =  xbd%HK2_T2 * HKREF
           xbd%VS2(4,3) =  xbd%HK2_D2 * HKREF
           xbd%VS2(4,4) =( xbd%HK2_U2 * HKREF  +  SENS/UEREF )
     &  *xbd%U2_UEI
           xbd%VSREZ(4) = -(HKREF**2)*(xbd%HK2 / HKREF - 1.0)
     &                     - SENS*(xbd%U2  / UEREF - 1.0)
C
          ENDIF
C
C-------- solve Newton system for current "2" station
          CALL GAUSS(4,4,xbd%VS2,xbd%VSREZ,1)
C
C-------- determine max changes and underrelax if necessary
C-------- (added Ue clamp   MD  3 Apr 03)
          DMAX = MAX( ABS(xbd%VSREZ(2)/THI),
     &                ABS(xbd%VSREZ(3)/DSI),
     &                ABS(xbd%VSREZ(4)/UEI)  )
          IF(IBL.GE.xfd%ITRAN(IS)) DMAX = MAX(DMAX,ABS(xbd%VSREZ(1)/(10
     &  .0*CTI)))
C
          xfd%RLX = 1.0
          IF(DMAX.GT.0.3) xfd%RLX = 0.3/DMAX
C
C-------- update as usual
          IF(IBL.LT.xfd%ITRAN(IS)) AMI = AMI + xfd%RLX*xbd%VSREZ(1)
          IF(IBL.GE.xfd%ITRAN(IS)) CTI = CTI + xfd%RLX*xbd%VSREZ(1)
          THI = THI + xfd%RLX*xbd%VSREZ(2)
          DSI = DSI + xfd%RLX*xbd%VSREZ(3)
          UEI = UEI + xfd%RLX*xbd%VSREZ(4)
C
C-------- eliminate absurd transients
          IF(IBL.GE.xfd%ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.xfd%IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*xbd%HSTINV / (xbd%GM1BL*(1.0 - 0.5*UEI*UEI
     &  *xbd%HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.DEPS) GO TO 110
C
  100   CONTINUE
C
        IF (.NOT. xfd%SILENT_MODE) WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHDU: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN
          IF(IBL.LE.xfd%IBLTE(IS)) THEN
           THI = xfd%THET(IBM,IS) * (xfd%XSSI(IBL,IS)/xfd%XSSI(IBM,IS))*
     &  *0.5
           DSI = xfd%DSTR(IBM,IS) * (xfd%XSSI(IBL,IS)/xfd%XSSI(IBM,IS))*
     &  *0.5
           UEI = xfd%UEDG(IBM,IS)
          ELSE IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
           UEI = xfd%UEDG(IBM,IS)
          ELSE
           THI = xfd%THET(IBM,IS)
           RATLEN = (xfd%XSSI(IBL,IS)-xfd%XSSI(IBM,IS)) / (10.0
     &  *xfd%DSTR(IBM,IS))
           DSI = (xfd%DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
           UEI = xfd%UEDG(IBM,IS)
          ENDIF
          IF(IBL.EQ.xfd%ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.xfd%ITRAN(IS)) CTI = xfd%CTAU(IBM,IS)
         ENDIF
C
 109     CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN(xbd)
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.xbd%SIMI) .AND. (.NOT.xbd%TURB)) THEN
          CALL TRCHEK(xfd,xbd)
C         DP mod: check for infinite loop condition
          IF (xfd%XFOIL_FAIL) RETURN
          AMI = xbd%AMPL2
          IF(     xbd%TRAN) xfd%ITRAN(IS) = IBL
          IF(.NOT.xbd%TRAN) xfd%ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.xfd%ITRAN(IS)) CALL BLVAR(bld,xbd,1)
         IF(IBL.GE.xfd%ITRAN(IS)) CALL BLVAR(bld,xbd,2)
         IF(xbd%WAKE) CALL BLVAR(bld,xbd,3)
C
         IF(IBL.LT.xfd%ITRAN(IS)) CALL BLMID(bld,xbd,1)
         IF(IBL.GE.xfd%ITRAN(IS)) CALL BLMID(bld,xbd,2)
         IF(xbd%WAKE) CALL BLMID(bld,xbd,3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
        SENS = SENNEW
C
C------ store primary variables
        IF(IBL.LT.xfd%ITRAN(IS)) xfd%CTAU(IBL,IS) = AMI
        IF(IBL.GE.xfd%ITRAN(IS)) xfd%CTAU(IBL,IS) = CTI
        xfd%THET(IBL,IS) = THI
        xfd%DSTR(IBL,IS) = DSI
        xfd%UEDG(IBL,IS) = UEI
        xfd%MASS(IBL,IS) = DSI*UEI
        xfd%TAU(IBL,IS)  = 0.5*xbd%R2*xbd%U2*xbd%U2*xbd%CF2
        xfd%DIS(IBL,IS)  =     xbd%R2*xbd%U2*xbd%U2*xbd%U2*xbd%DI2
     &  *xbd%HS2*0.5
        xfd%CTQ(IBL,IS)  = xbd%CQ2
        xfd%DELT(IBL,IS) = xbd%DE2
        xfd%TSTR(IBL,IS) = xbd%HS2*xbd%T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN(xbd)
C       DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
        call com2_to_com1(xbd)
C        DO 310 ICOM=1, NCOM
C          COM1(ICOM) = COM2(ICOM)
C  310   CONTINUE
C
C
C------ turbulent intervals will follow transition interval or TE
        IF(xbd%TRAN .OR. IBL.EQ.xfd%IBLTE(IS)) THEN
         xbd%TURB = .TRUE.
C
C------- save transition location
         xfd%TFORCE(IS) = xbd%TRFORC
         xfd%XSSITR(IS) = xbd%XT
        ENDIF
C
        xbd%TRAN = .FALSE.
C
 1000 CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets Ue from inviscid Ue plus all source influence
C
C===================================================================70
      SUBROUTINE UESET(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      DO 1 IS=1, 2
        DO 10 IBL=2, xfd%NBL(IS)
          I = xfd%IPAN(IBL,IS)
C
          DUI = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, xfd%NBL(JS)
              J  = xfd%IPAN(JBL,JS)
              UE_M = -xfd%VTI(IBL,IS)*xfd%VTI(JBL,JS)*xfd%DIJ(I,J)
              DUI = DUI + UE_M*xfd%MASS(JBL,JS)
 1000       CONTINUE
  100     CONTINUE
C
          xfd%UEDG(IBL,IS) = xfd%UINV(IBL,IS) + DUI
C
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets up the BL Newton system coefficients
C     for the current BL variables and the edge
C     velocities received from SETUP. The local
C     BL system coefficients are then
C     incorporated into the global Newton system.  
C
C===================================================================70
      SUBROUTINE SETBL(xfd,bld,xbd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd

      REAL*8 USAV(IVX,2)
      REAL*8 :: U1_M(2*IVX), U2_M(2*IVX)
      REAL*8 :: D1_M(2*IVX), D2_M(2*IVX)
      REAL*8 :: ULE1_M(2*IVX), ULE2_M(2*IVX)
      REAL*8 :: UTE1_M(2*IVX), UTE2_M(2*IVX)
      REAL*8 :: MA_CLMR, MSQ_CLMR, MDI

C     DP mod: set explicitly (otherwise initialized as 0 here)
      HVRAT = 0.0
C
C---- set the CL used to define Mach, Reynolds numbers
      IF(xfd%LALFA) THEN
       CLMR = xfd%CL
      ELSE
       CLMR = xfd%CLSPEC
      ENDIF
C
C---- set current MINF(CL)
      CALL MRCL(xfd,CLMR,MA_CLMR,RE_CLMR)
      MSQ_CLMR = 2.0*xfd%MINF*MA_CLMR
C
C---- set compressibility parameter TKLAM and derivative TK_MSQ
      CALL COMSET(xfd)
C
C---- set gas constant (= Cp/Cv)
      xbd%GAMBL = xfd%GAMMA
      xbd%GM1BL = xfd%GAMM1
C
C---- set parameters for compressibility correction
      xbd%QINFBL = xfd%QINF
      xbd%TKBL    = xfd%TKLAM
      xbd%TKBL_MS = xfd%TKL_MSQ
C
C---- stagnation density and 1/enthalpy
      xbd%RSTBL    = (1.0 + 0.5*xbd%GM1BL*xfd%MINF**2) ** (1.0/xbd%GM1BL
     &  )
      xbd%RSTBL_MS = 0.5*xbd%RSTBL/(1.0 + 0.5*xbd%GM1BL*xfd%MINF**2)
C
      xbd%HSTINV    = xbd%GM1BL*(xfd%MINF/xbd%QINFBL)**2 / (1.0 + 0.5
     &  *xbd%GM1BL*xfd%MINF**2)
      xbd%HSTINV_MS = xbd%GM1BL*( 1.0/xbd%QINFBL)**2 / (1.0 + 0.5
     &  *xbd%GM1BL*xfd%MINF**2)
     &                - 0.5*xbd%GM1BL*xbd%HSTINV / (1.0 + 0.5*xbd%GM1BL
     &  *xfd%MINF**2)
C
C---- set Reynolds number based on freestream density, velocity, viscosity
      HERAT    = 1.0 - 0.5*xbd%QINFBL**2*xbd%HSTINV
      HERAT_MS =     - 0.5*xbd%QINFBL**2*xbd%HSTINV_MS
C
      xbd%REYBL    = xfd%REINF * SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT
     &  +HVRAT)
      xbd%REYBL_RE =         SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
      xbd%REYBL_MS = xbd%REYBL * (1.5/HERAT - 1.0/(HERAT+HVRAT))
     &  *HERAT_MS
C
!FIXME: xfd%IDAMP always set to 0 here.  In xfoil, can be set to 1 by user.
      xbd%IDAMPV = xfd%IDAMP
C
C---- save TE thickness
      xbd%DWTE = xfd%WGAP(1)
C
      IF(.NOT.xfd%LBLINI) THEN
C----- initialize BL by marching with Ue (fudge at separation)
C      DP mod: added SILENT_MODE option
       IF (.NOT. xfd%SILENT_MODE) THEN
         WRITE(*,*)
         WRITE(*,*) 'Initializing BL ...'
       ENDIF
       CALL MRCHUE(xfd,bld,xbd)
C      DP mod: check for infinite loop condition
       IF (xfd%XFOIL_FAIL) RETURN
       xfd%LBLINI = .TRUE.
      ENDIF
C
      IF (.NOT. xfd%SILENT_MODE) WRITE(*,*)
C
C---- march BL with current Ue and Ds to establish transition
      CALL MRCHDU(xfd,bld,xbd)
C     DP mod: check for infinite loop condition
      IF (xfd%XFOIL_FAIL) RETURN
C
      DO 5 IS=1, 2
        DO 6 IBL=2, xfd%NBL(IS)
          USAV(IBL,IS) = xfd%UEDG(IBL,IS)
    6   CONTINUE
    5 CONTINUE
C
      CALL UESET(xfd)
C
      DO 7 IS=1, 2
        DO 8 IBL=2, xfd%NBL(IS)
          TEMP = USAV(IBL,IS)
          USAV(IBL,IS) = xfd%UEDG(IBL,IS)
          xfd%UEDG(IBL,IS) = TEMP
    8   CONTINUE
    7 CONTINUE
C
      ILE1 = xfd%IPAN(2,1)
      ILE2 = xfd%IPAN(2,2)
      ITE1 = xfd%IPAN(xfd%IBLTE(1),1)
      ITE2 = xfd%IPAN(xfd%IBLTE(2),2)
C
      JVTE1 = xfd%ISYS(xfd%IBLTE(1),1)
      JVTE2 = xfd%ISYS(xfd%IBLTE(2),2)
C
      DULE1 = xfd%UEDG(2,1) - USAV(2,1)
      DULE2 = xfd%UEDG(2,2) - USAV(2,2)
C
C---- set LE and TE Ue sensitivities wrt all m values
      DO 10 JS=1, 2
        DO 110 JBL=2, xfd%NBL(JS)
          J  = xfd%IPAN(JBL,JS)
          JV = xfd%ISYS(JBL,JS)
          ULE1_M(JV) = -xfd%VTI(       2,1)*xfd%VTI(JBL,JS)*xfd%DIJ(ILE1
     &  ,J)
          ULE2_M(JV) = -xfd%VTI(       2,2)*xfd%VTI(JBL,JS)*xfd%DIJ(ILE2
     &  ,J)
          UTE1_M(JV) = -xfd%VTI(xfd%IBLTE(1),1)*xfd%VTI(JBL,JS)
     &  *xfd%DIJ(ITE1,J)
          UTE2_M(JV) = -xfd%VTI(xfd%IBLTE(2),2)*xfd%VTI(JBL,JS)
     &  *xfd%DIJ(ITE2,J)
  110   CONTINUE
   10 CONTINUE
C
      ULE1_A = xfd%UINV_A(2,1)
      ULE2_A = xfd%UINV_A(2,2)

C      DP mod: TINDEX not used (it is for storing polars in Xfoil)
C      TINDEX(1) = 0.0
C      TINDEX(2) = 0.0
C
C**** Go over each boundary layer/wake
      DO 2000 IS=1, 2
C
C---- there is no station "1" at similarity, so zero everything out
      DO 20 JS=1, 2
        DO 210 JBL=2, xfd%NBL(JS)
          JV = xfd%ISYS(JBL,JS)
          U1_M(JV) = 0.
          D1_M(JV) = 0.
  210   CONTINUE
   20 CONTINUE
      U1_A = 0.
      D1_A = 0.
C
      DUE1 = 0.
      DDS1 = 0.
C
C---- similarity station pressure gradient parameter  x/u du/dx
      IBL = 2
      xbd%BULE = 1.0
C
      xbd%AMCRIT = xfd%ACRIT
C
C---- set forced transition arc length position
      CALL XIFSET(xfd,xbd,IS)
C
      xbd%TRAN = .FALSE.
      xbd%TURB = .FALSE.
C
C**** Sweep downstream setting up BL equation linearizations
      DO 1000 IBL=2, xfd%NBL(IS)
C
      IV  = xfd%ISYS(IBL,IS)
C
      xbd%SIMI = IBL.EQ.2
      xbd%WAKE = IBL.GT.xfd%IBLTE(IS)
      xbd%TRAN = IBL.EQ.xfd%ITRAN(IS)
      xbd%TURB = IBL.GT.xfd%ITRAN(IS)
C
      I = xfd%IPAN(IBL,IS)
C
C---- set primary variables for current station
      XSI = xfd%XSSI(IBL,IS)
      IF(IBL.LT.xfd%ITRAN(IS)) AMI = xfd%CTAU(IBL,IS)
      IF(IBL.GE.xfd%ITRAN(IS)) CTI = xfd%CTAU(IBL,IS)
      UEI = xfd%UEDG(IBL,IS)
      THI = xfd%THET(IBL,IS)
      MDI = xfd%MASS(IBL,IS)
C
      DSI = MDI/UEI
C
      IF(xbd%WAKE) THEN
       IW = IBL - xfd%IBLTE(IS)
       DSWAKI = xfd%WGAP(IW)
      ELSE
       DSWAKI = 0.
      ENDIF
C
C---- set derivatives of DSI (= D2)
      D2_M2 =  1.0/UEI
      D2_U2 = -DSI/UEI
C
      DO 30 JS=1, 2
        DO 310 JBL=2, xfd%NBL(JS)
          J  = xfd%IPAN(JBL,JS)
          JV = xfd%ISYS(JBL,JS)
          U2_M(JV) = -xfd%VTI(IBL,IS)*xfd%VTI(JBL,JS)*xfd%DIJ(I,J)
          D2_M(JV) = D2_U2*U2_M(JV)
  310   CONTINUE
   30 CONTINUE
      D2_M(IV) = D2_M(IV) + D2_M2
C
      U2_A = xfd%UINV_A(IBL,IS)
      D2_A = D2_U2*U2_A
C
C---- "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
      DUE2 = xfd%UEDG(IBL,IS) - USAV(IBL,IS)
      DDS2 = D2_U2*DUE2
C
      CALL BLPRV(xbd,XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
      CALL BLKIN(xbd)
C
C---- check for transition and set TRAN, XT, etc. if found
      IF(xbd%TRAN) THEN
        CALL TRCHEK(xfd,xbd)
C       DP mod: check for infinite loop condition
        IF (xfd%XFOIL_FAIL) RETURN
        AMI = xbd%AMPL2
      ENDIF
C     DP mod: added SILENT_MODE option
      IF(IBL.EQ.xfd%ITRAN(IS) .AND. .NOT.xbd%TRAN .AND. .NOT
     &  .xfd%SILENT_MODE) THEN
       WRITE(*,*) 'SETBL: Xtr???  n1 n2: ', xbd%AMPL1, xbd%AMPL2
      ENDIF
C
C---- assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
C     at the previous "1" station and the current "2" station
C
      IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
C
C----- define quantities at start of wake, adding TE base thickness to Dstar
       TTE = xfd%THET(xfd%IBLTE(1),1) + xfd%THET(xfd%IBLTE(2),2)
       DTE = xfd%DSTR(xfd%IBLTE(1),1) + xfd%DSTR(xfd%IBLTE(2),2) +
     &   xfd%ANTE
       CTE = ( xfd%CTAU(xfd%IBLTE(1),1)*xfd%THET(xfd%IBLTE(1),1)
     &       + xfd%CTAU(xfd%IBLTE(2),2)*xfd%THET(xfd%IBLTE(2),2) ) / TTE
     &  
       CALL TESYS(bld,xbd,CTE,TTE,DTE)
C
       TTE_TTE1 = 1.0
       TTE_TTE2 = 1.0
       DTE_MTE1 =               1.0 / xfd%UEDG(xfd%IBLTE(1),1)
       DTE_UTE1 = -xfd%DSTR(xfd%IBLTE(1),1) / xfd%UEDG(xfd%IBLTE(1),1)
       DTE_MTE2 =               1.0 / xfd%UEDG(xfd%IBLTE(2),2)
       DTE_UTE2 = -xfd%DSTR(xfd%IBLTE(2),2) / xfd%UEDG(xfd%IBLTE(2),2)
       CTE_CTE1 = xfd%THET(xfd%IBLTE(1),1)/TTE
       CTE_CTE2 = xfd%THET(xfd%IBLTE(2),2)/TTE
       CTE_TTE1 = (xfd%CTAU(xfd%IBLTE(1),1) - CTE)/TTE
       CTE_TTE2 = (xfd%CTAU(xfd%IBLTE(2),2) - CTE)/TTE
C
C----- re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
       DO 35 JS=1, 2
         DO 350 JBL=2, xfd%NBL(JS)
           J  = xfd%IPAN(JBL,JS)
           JV = xfd%ISYS(JBL,JS)
           D1_M(JV) = DTE_UTE1*UTE1_M(JV) + DTE_UTE2*UTE2_M(JV)
  350    CONTINUE
   35  CONTINUE
       D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
       D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
C
C----- "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
       DUE1 = 0.
       DDS1 = DTE_UTE1*(xfd%UEDG(xfd%IBLTE(1),1) - USAV(xfd%IBLTE(1),1))
     &  
     &      + DTE_UTE2*(xfd%UEDG(xfd%IBLTE(2),2) - USAV(xfd%IBLTE(2),2))
     &  
C
      ELSE
C
       CALL BLSYS(bld,xbd)
C
      ENDIF
C
C
C---- Save wall shear and equil. max shear coefficient for plotting output
      xfd%TAU(IBL,IS) = 0.5*xbd%R2*xbd%U2*xbd%U2*xbd%CF2
      xfd%DIS(IBL,IS) =     xbd%R2*xbd%U2*xbd%U2*xbd%U2*xbd%DI2*xbd%HS2
     &  *0.5
      xfd%CTQ(IBL,IS) = xbd%CQ2
      xfd%DELT(IBL,IS) = xbd%DE2
      xfd%USLP(IBL,IS) = 1.60/(1.0+xbd%US2)
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c      IF(WAKE) THEN
c        ALD = DLCON
c      ELSE
c       ALD = 1.0
c      ENDIF
cC
c      IF(TURB .AND. .NOT.WAKE) THEN
c        GCC = GCCON
c        HKC     = HK2 - 1.0 - GCC/RT2
c        IF(HKC .LT. 0.01) THEN
c         HKC = 0.01
c        ENDIF
c       ELSE
c        HKC = HK2 - 1.0
c       ENDIF
cC
c       HR = HKC     / (GACON*ALD*HK2)
c       UQ = (0.5*CF2 - HR**2) / (GBCON*D2)
cC
c       IF(TURB) THEN
c        IBLP = MIN(IBL+1,NBL(IS))
c        IBLM = MAX(IBL-1,2      )
c        DXSSI = XSSI(IBLP,IS) - XSSI(IBLM,IS)
c        IF(DXXSI.EQ.0.0) DXSSI = 1.0
c        GUXD(IBL,IS) = -LOG(UEDG(IBLP,IS)/UEDG(IBLM,IS)) / DXSSI
c        GUXQ(IBL,IS) = -UQ
c       ELSE
c        GUXD(IBL,IS) = 0.0
c        GUXQ(IBL,IS) = 0.0
c       ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C---- set XI sensitivities wrt LE Ue changes
      IF(IS.EQ.1) THEN
       XI_ULE1 =  xfd%SST_GO
       XI_ULE2 = -xfd%SST_GP
      ELSE
       XI_ULE1 = -xfd%SST_GO
       XI_ULE2 =  xfd%SST_GP
      ENDIF
C
C---- stuff BL system coefficients into main Jacobian matrix
C
      DO 40 JV=1, xfd%NSYS
        xfd%VM(1,JV,IV) = xbd%VS1(1,3)*D1_M(JV) + xbd%VS1(1,4)*U1_M(JV)
     &              + xbd%VS2(1,3)*D2_M(JV) + xbd%VS2(1,4)*U2_M(JV)
     &              + (xbd%VS1(1,5) + xbd%VS2(1,5) + xbd%VSX(1))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   40 CONTINUE
C
      xfd%VB(1,1,IV) = xbd%VS1(1,1)
      xfd%VB(1,2,IV) = xbd%VS1(1,2)
C
      xfd%VA(1,1,IV) = xbd%VS2(1,1)
      xfd%VA(1,2,IV) = xbd%VS2(1,2)
C
      IF(xfd%LALFA) THEN
       xfd%VDEL(1,2,IV) = xbd%VSR(1)*RE_CLMR + xbd%VSM(1)*MSQ_CLMR
      ELSE
       xfd%VDEL(1,2,IV) = 
     &       (xbd%VS1(1,4)*U1_A + xbd%VS1(1,3)*D1_A)
     &     + (xbd%VS2(1,4)*U2_A + xbd%VS2(1,3)*D2_A)
     &     + (xbd%VS1(1,5) + xbd%VS2(1,5) + xbd%VSX(1))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      xfd%VDEL(1,1,IV) = xbd%VSREZ(1)
     &   + (xbd%VS1(1,4)*DUE1 + xbd%VS1(1,3)*DDS1)
     &   + (xbd%VS2(1,4)*DUE2 + xbd%VS2(1,3)*DDS2)
     &   + (xbd%VS1(1,5) + xbd%VS2(1,5) + xbd%VSX(1))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 50 JV=1, xfd%NSYS
        xfd%VM(2,JV,IV) = xbd%VS1(2,3)*D1_M(JV) + xbd%VS1(2,4)*U1_M(JV)
     &              + xbd%VS2(2,3)*D2_M(JV) + xbd%VS2(2,4)*U2_M(JV)
     &              + (xbd%VS1(2,5) + xbd%VS2(2,5) + xbd%VSX(2))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   50 CONTINUE
C
      xfd%VB(2,1,IV)  = xbd%VS1(2,1)
      xfd%VB(2,2,IV)  = xbd%VS1(2,2)
C
      xfd%VA(2,1,IV) = xbd%VS2(2,1)
      xfd%VA(2,2,IV) = xbd%VS2(2,2)
C
      IF(xfd%LALFA) THEN
       xfd%VDEL(2,2,IV) = xbd%VSR(2)*RE_CLMR + xbd%VSM(2)*MSQ_CLMR
      ELSE
       xfd%VDEL(2,2,IV) = 
     &       (xbd%VS1(2,4)*U1_A + xbd%VS1(2,3)*D1_A)
     &     + (xbd%VS2(2,4)*U2_A + xbd%VS2(2,3)*D2_A)
     &     + (xbd%VS1(2,5) + xbd%VS2(2,5) + xbd%VSX(2))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      xfd%VDEL(2,1,IV) = xbd%VSREZ(2)
     &   + (xbd%VS1(2,4)*DUE1 + xbd%VS1(2,3)*DDS1)
     &   + (xbd%VS2(2,4)*DUE2 + xbd%VS2(2,3)*DDS2)
     &   + (xbd%VS1(2,5) + xbd%VS2(2,5) + xbd%VSX(2))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 60 JV=1, xfd%NSYS
        xfd%VM(3,JV,IV) = xbd%VS1(3,3)*D1_M(JV) + xbd%VS1(3,4)*U1_M(JV)
     &              + xbd%VS2(3,3)*D2_M(JV) + xbd%VS2(3,4)*U2_M(JV)
     &              + (xbd%VS1(3,5) + xbd%VS2(3,5) + xbd%VSX(3))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   60 CONTINUE
C
      xfd%VB(3,1,IV) = xbd%VS1(3,1)
      xfd%VB(3,2,IV) = xbd%VS1(3,2)
C
      xfd%VA(3,1,IV) = xbd%VS2(3,1)
      xfd%VA(3,2,IV) = xbd%VS2(3,2)
C
      IF(xfd%LALFA) THEN
       xfd%VDEL(3,2,IV) = xbd%VSR(3)*RE_CLMR + xbd%VSM(3)*MSQ_CLMR
      ELSE
       xfd%VDEL(3,2,IV) = 
     &       (xbd%VS1(3,4)*U1_A + xbd%VS1(3,3)*D1_A)
     &     + (xbd%VS2(3,4)*U2_A + xbd%VS2(3,3)*D2_A)
     &     + (xbd%VS1(3,5) + xbd%VS2(3,5) + xbd%VSX(3))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      xfd%VDEL(3,1,IV) = xbd%VSREZ(3)
     &   + (xbd%VS1(3,4)*DUE1 + xbd%VS1(3,3)*DDS1)
     &   + (xbd%VS2(3,4)*DUE2 + xbd%VS2(3,3)*DDS2)
     &   + (xbd%VS1(3,5) + xbd%VS2(3,5) + xbd%VSX(3))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      IF(IBL.EQ.xfd%IBLTE(IS)+1) THEN
C
C----- redefine coefficients for TTE, DTE, etc
       xfd%VZ(1,1)    = xbd%VS1(1,1)*CTE_CTE1
       xfd%VZ(1,2)    = xbd%VS1(1,1)*CTE_TTE1 + xbd%VS1(1,2)*TTE_TTE1
       xfd%VB(1,1,IV) = xbd%VS1(1,1)*CTE_CTE2
       xfd%VB(1,2,IV) = xbd%VS1(1,1)*CTE_TTE2 + xbd%VS1(1,2)*TTE_TTE2
C
       xfd%VZ(2,1)    = xbd%VS1(2,1)*CTE_CTE1
       xfd%VZ(2,2)    = xbd%VS1(2,1)*CTE_TTE1 + xbd%VS1(2,2)*TTE_TTE1
       xfd%VB(2,1,IV) = xbd%VS1(2,1)*CTE_CTE2
       xfd%VB(2,2,IV) = xbd%VS1(2,1)*CTE_TTE2 + xbd%VS1(2,2)*TTE_TTE2
C
       xfd%VZ(3,1)    = xbd%VS1(3,1)*CTE_CTE1
       xfd%VZ(3,2)    = xbd%VS1(3,1)*CTE_TTE1 + xbd%VS1(3,2)*TTE_TTE1
       xfd%VB(3,1,IV) = xbd%VS1(3,1)*CTE_CTE2
       xfd%VB(3,2,IV) = xbd%VS1(3,1)*CTE_TTE2 + xbd%VS1(3,2)*TTE_TTE2
C
      ENDIF
C
C---- turbulent intervals will follow if currently at transition interval
      IF(xbd%TRAN) THEN
        xbd%TURB = .TRUE.
C
C------ save transition location
        xfd%ITRAN(IS) = IBL
        xfd%TFORCE(IS) = xbd%TRFORC
        xfd%XSSITR(IS) = xbd%XT
C
C------ interpolate airfoil geometry to find transition x/c
C-      (for user output)
        IF(IS.EQ.1) THEN
         STR = xfd%SST - xbd%XT
        ELSE
         STR = xfd%SST + xbd%XT
        ENDIF
        CHX = xfd%XTE - xfd%XLE
        CHY = xfd%YTE - xfd%YLE
        CHSQ = CHX**2 + CHY**2
        XTR = SEVAL(STR,xfd%X,xfd%XP,xfd%S,xfd%N)
        YTR = SEVAL(STR,xfd%Y,xfd%YP,xfd%S,xfd%N)
        xfd%XOCTR(IS) = ((XTR-xfd%XLE)*CHX + (YTR-xfd%YLE)*CHY)/CHSQ
        xfd%YOCTR(IS) = ((YTR-xfd%YLE)*CHX - (XTR-xfd%XLE)*CHY)/CHSQ
      ENDIF
C
      xbd%TRAN = .FALSE.
C
      IF(IBL.EQ.xfd%IBLTE(IS)) THEN
C----- set "2" variables at TE to wake correlations for next station
C
       xbd%TURB = .TRUE.
       xbd%WAKE = .TRUE.
       CALL BLVAR(bld,xbd,3)
       CALL BLMID(bld,xbd,3)
      ENDIF
C
      DO 80 JS=1, 2
        DO 810 JBL=2, xfd%NBL(JS)
          JV = xfd%ISYS(JBL,JS)
          U1_M(JV) = U2_M(JV)
          D1_M(JV) = D2_M(JV)
  810   CONTINUE
   80 CONTINUE
C
      U1_A = U2_A
      D1_A = D2_A
C
      DUE1 = DUE2
      DDS1 = DDS2
C
C      DP mod: TINDEX not needed (used for storing polars in Xfoil)
C      IF(IBL .EQ. ITRAN(IS) .AND. X2 .GT. X1) THEN
C       IF(IS.EQ.1) THEN
C        TINDEX(IS) = FLOAT(IST-ITRAN(IS)+3) - (XT-X1)/(X2-X1)
C       ELSE
C        TINDEX(IS) = FLOAT(IST+ITRAN(IS)-2) + (XT-X1)/(X2-X1)
C       ENDIF
C      ENDIF
C      
C---- set BL variables for next station
C     DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
      call com2_to_com1(xbd)
C      DO 190 ICOM=1, NCOM
C        COM1(ICOM) = COM2(ICOM)
C  190 CONTINUE
C
C---- next streamwise station
 1000 CONTINUE
C 
C     DP mod: added SILENT_MODE option
C     DP mod: change write precision
      IF (.NOT. xfd%SILENT_MODE) THEN
        IF(xfd%TFORCE(IS)) THEN
         WRITE(*,9100) IS,xfd%XOCTR(IS),xfd%ITRAN(IS)
 9100    FORMAT(1X,'Side',I2,' forced transition at x/c = ',F15.12,I5)
        ELSE
         WRITE(*,9200) IS,xfd%XOCTR(IS),xfd%ITRAN(IS)
 9200    FORMAT(1X,'Side',I2,'  free  transition at x/c = ',F15.12,I5)
        ENDIF
      ENDIF
C
C---- next airfoil side
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Custom solver for coupled viscous-inviscid Newton system:
C
C        A  |  |  .  |  |  .  |    d       R       S
C        B  A  |  .  |  |  .  |    d       R       S
C        |  B  A  .  |  |  .  |    d       R       S
C        .  .  .  .  |  |  .  |    d   =   R - dRe S
C        |  |  |  B  A  |  .  |    d       R       S
C        |  Z  |  |  B  A  .  |    d       R       S
C        .  .  .  .  .  .  .  |    d       R       S
C        |  |  |  |  |  |  B  A    d       R       S
C
C       A, B, Z  3x3  blocks containing linearized BL equation coefficients
C       |        3x1  vectors containing mass defect influence 
C                     coefficients on Ue
C       d        3x1  unknown vectors (Newton deltas for Ctau, Theta, m)
C       R        3x1  residual vectors
C       S        3x1  Re influence vectors
C
C===================================================================70
      SUBROUTINE BLSOLV(xfd)

      use my_equivalence, only : my_equiv_3_2, my_equiv_3_1
      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      IVTE1 = xfd%ISYS(xfd%IBLTE(1),1)
C
      VACC1 = xfd%VACCEL
      VACC2 = xfd%VACCEL * 2.0 / (xfd%S(xfd%N) - xfd%S(1))
      VACC3 = xfd%VACCEL * 2.0 / (xfd%S(xfd%N) - xfd%S(1))
C
      DO 1000 IV=1, xfd%NSYS
C
        IVP = IV + 1
C
C====== Invert VA(IV) block
C
C------ normalize first row
        PIVOT = 1.0 / xfd%VA(1,1,IV)
        xfd%VA(1,2,IV) = xfd%VA(1,2,IV) * PIVOT
        DO 10 L=IV, xfd%NSYS
          xfd%VM(1,L,IV) = xfd%VM(1,L,IV)*PIVOT
   10   CONTINUE
        xfd%VDEL(1,1,IV) = xfd%VDEL(1,1,IV)*PIVOT
        xfd%VDEL(1,2,IV) = xfd%VDEL(1,2,IV)*PIVOT
C
C------ eliminate lower first column in VA block
        DO 15 K=2, 3
          VTMP = xfd%VA(K,1,IV)
          xfd%VA(K,2,IV) = xfd%VA(K,2,IV) - VTMP*xfd%VA(1,2,IV)
          DO 150 L=IV, xfd%NSYS
            xfd%VM(K,L,IV) = xfd%VM(K,L,IV) - VTMP*xfd%VM(1,L,IV)
  150     CONTINUE
          xfd%VDEL(K,1,IV) = xfd%VDEL(K,1,IV) - VTMP*xfd%VDEL(1,1,IV)
          xfd%VDEL(K,2,IV) = xfd%VDEL(K,2,IV) - VTMP*xfd%VDEL(1,2,IV)
   15   CONTINUE
C
C
C------ normalize second row
        PIVOT = 1.0 / xfd%VA(2,2,IV)
        DO 20 L=IV, xfd%NSYS
          xfd%VM(2,L,IV) = xfd%VM(2,L,IV)*PIVOT
   20   CONTINUE
        xfd%VDEL(2,1,IV) = xfd%VDEL(2,1,IV)*PIVOT
        xfd%VDEL(2,2,IV) = xfd%VDEL(2,2,IV)*PIVOT
C
C------ eliminate lower second column in VA block
        K = 3
        VTMP = xfd%VA(K,2,IV)
        DO 250 L=IV, xfd%NSYS
          xfd%VM(K,L,IV) = xfd%VM(K,L,IV) - VTMP*xfd%VM(2,L,IV)
  250   CONTINUE
        xfd%VDEL(K,1,IV) = xfd%VDEL(K,1,IV) - VTMP*xfd%VDEL(2,1,IV)
        xfd%VDEL(K,2,IV) = xfd%VDEL(K,2,IV) - VTMP*xfd%VDEL(2,2,IV)
C
C
C------ normalize third row
        PIVOT = 1.0/xfd%VM(3,IV,IV)
        DO 350 L=IVP, xfd%NSYS
          xfd%VM(3,L,IV) = xfd%VM(3,L,IV)*PIVOT
  350   CONTINUE
        xfd%VDEL(3,1,IV) = xfd%VDEL(3,1,IV)*PIVOT
        xfd%VDEL(3,2,IV) = xfd%VDEL(3,2,IV)*PIVOT
C
C
C------ eliminate upper third column in VA block
        VTMP1 = xfd%VM(1,IV,IV)
        VTMP2 = xfd%VM(2,IV,IV)
        DO 450 L=IVP, xfd%NSYS
          xfd%VM(1,L,IV) = xfd%VM(1,L,IV) - VTMP1*xfd%VM(3,L,IV)
          xfd%VM(2,L,IV) = xfd%VM(2,L,IV) - VTMP2*xfd%VM(3,L,IV)
  450   CONTINUE
        xfd%VDEL(1,1,IV) = xfd%VDEL(1,1,IV) - VTMP1*xfd%VDEL(3,1,IV)
        xfd%VDEL(2,1,IV) = xfd%VDEL(2,1,IV) - VTMP2*xfd%VDEL(3,1,IV)
        xfd%VDEL(1,2,IV) = xfd%VDEL(1,2,IV) - VTMP1*xfd%VDEL(3,2,IV)
        xfd%VDEL(2,2,IV) = xfd%VDEL(2,2,IV) - VTMP2*xfd%VDEL(3,2,IV)
C
C------ eliminate upper second column in VA block
        VTMP = xfd%VA(1,2,IV)
        DO 460 L=IVP, xfd%NSYS
          xfd%VM(1,L,IV) = xfd%VM(1,L,IV) - VTMP*xfd%VM(2,L,IV)
  460   CONTINUE
        xfd%VDEL(1,1,IV) = xfd%VDEL(1,1,IV) - VTMP*xfd%VDEL(2,1,IV)
        xfd%VDEL(1,2,IV) = xfd%VDEL(1,2,IV) - VTMP*xfd%VDEL(2,2,IV)
C
C
        IF(IV.EQ.xfd%NSYS) GO TO 1000
C
C====== Eliminate VB(IV+1) block, rows  1 -> 3
        DO 50 K=1, 3
          VTMP1 = xfd%VB(K, 1,IVP)
          VTMP2 = xfd%VB(K, 2,IVP)
          VTMP3 = xfd%VM(K,IV,IVP)
          DO 510 L=IVP, xfd%NSYS
            xfd%VM(K,L,IVP) = xfd%VM(K,L,IVP)
     &        - (  VTMP1*xfd%VM(1,L,IV)
     &           + VTMP2*xfd%VM(2,L,IV)
     &           + VTMP3*xfd%VM(3,L,IV) )
  510     CONTINUE
          xfd%VDEL(K,1,IVP) = xfd%VDEL(K,1,IVP)
     &        - (  VTMP1*xfd%VDEL(1,1,IV)
     &           + VTMP2*xfd%VDEL(2,1,IV)
     &           + VTMP3*xfd%VDEL(3,1,IV) )
          xfd%VDEL(K,2,IVP) = xfd%VDEL(K,2,IVP)
     &        - (  VTMP1*xfd%VDEL(1,2,IV)
     &           + VTMP2*xfd%VDEL(2,2,IV)
     &           + VTMP3*xfd%VDEL(3,2,IV) )
   50   CONTINUE
C
        IF(IV.EQ.IVTE1) THEN
C------- eliminate VZ block
         IVZ = xfd%ISYS(xfd%IBLTE(2)+1,2)
C
         DO 55 K=1, 3
           VTMP1 = xfd%VZ(K,1)
           VTMP2 = xfd%VZ(K,2)
           DO 515 L=IVP, xfd%NSYS
             xfd%VM(K,L,IVZ) = xfd%VM(K,L,IVZ)
     &         - (  VTMP1*xfd%VM(1,L,IV)
     &            + VTMP2*xfd%VM(2,L,IV) )
  515      CONTINUE
           xfd%VDEL(K,1,IVZ) = xfd%VDEL(K,1,IVZ)
     &         - (  VTMP1*xfd%VDEL(1,1,IV)
     &            + VTMP2*xfd%VDEL(2,1,IV) )
           xfd%VDEL(K,2,IVZ) = xfd%VDEL(K,2,IVZ)
     &         - (  VTMP1*xfd%VDEL(1,2,IV)
     &            + VTMP2*xfd%VDEL(2,2,IV) )
   55    CONTINUE
        ENDIF
C
        IF(IVP.EQ.xfd%NSYS) GO TO 1000
C
C====== Eliminate lower VM column
        DO 60 KV=IV+2, xfd%NSYS
          VTMP1 = xfd%VM(1,IV,KV)
          VTMP2 = xfd%VM(2,IV,KV)
          VTMP3 = xfd%VM(3,IV,KV)
C
          IF(ABS(VTMP1).GT.VACC1) THEN
          DO 610 L=IVP, xfd%NSYS
            xfd%VM(1,L,KV) = xfd%VM(1,L,KV) - VTMP1*xfd%VM(3,L,IV)
  610     CONTINUE
          xfd%VDEL(1,1,KV) = xfd%VDEL(1,1,KV) - VTMP1*xfd%VDEL(3,1,IV)
          xfd%VDEL(1,2,KV) = xfd%VDEL(1,2,KV) - VTMP1*xfd%VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP2).GT.VACC2) THEN
          DO 620 L=IVP, xfd%NSYS
            xfd%VM(2,L,KV) = xfd%VM(2,L,KV) - VTMP2*xfd%VM(3,L,IV)
  620     CONTINUE
          xfd%VDEL(2,1,KV) = xfd%VDEL(2,1,KV) - VTMP2*xfd%VDEL(3,1,IV)
          xfd%VDEL(2,2,KV) = xfd%VDEL(2,2,KV) - VTMP2*xfd%VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP3).GT.VACC3) THEN
          DO 630 L=IVP, xfd%NSYS
            xfd%VM(3,L,KV) = xfd%VM(3,L,KV) - VTMP3*xfd%VM(3,L,IV)
  630     CONTINUE
          xfd%VDEL(3,1,KV) = xfd%VDEL(3,1,KV) - VTMP3*xfd%VDEL(3,1,IV)
          xfd%VDEL(3,2,KV) = xfd%VDEL(3,2,KV) - VTMP3*xfd%VDEL(3,2,IV)
          ENDIF
C
   60   CONTINUE
C
 1000 CONTINUE
C
C
C
      DO 2000 IV=xfd%NSYS, 2, -1
C
C------ eliminate upper VM columns
        VTMP = xfd%VDEL(3,1,IV)
        DO 81 KV=IV-1, 1, -1
          xfd%VDEL(1,1,KV) = xfd%VDEL(1,1,KV) - xfd%VM(1,IV,KV)*VTMP
          xfd%VDEL(2,1,KV) = xfd%VDEL(2,1,KV) - xfd%VM(2,IV,KV)*VTMP
          xfd%VDEL(3,1,KV) = xfd%VDEL(3,1,KV) - xfd%VM(3,IV,KV)*VTMP
   81   CONTINUE
C
        VTMP = xfd%VDEL(3,2,IV)
        DO 82 KV=IV-1, 1, -1
          xfd%VDEL(1,2,KV) = xfd%VDEL(1,2,KV) - xfd%VM(1,IV,KV)*VTMP
          xfd%VDEL(2,2,KV) = xfd%VDEL(2,2,KV) - xfd%VM(2,IV,KV)*VTMP
          xfd%VDEL(3,2,KV) = xfd%VDEL(3,2,KV) - xfd%VM(3,IV,KV)*VTMP
   82   CONTINUE
C
 2000 CONTINUE
C
C     DP mod: copy from VA to UNEW to remove need for EQUIVALENCE
C       and from VA to U_AC to remove need for EQUIVALENCe
      DO 83 IV = 1, IVX
        call my_equiv_3_2(xfd%VA, xfd%UNEW, (/ 1, 1, 1 /), (/ 1, 1 /),
     &                   (/ IV, 1 /), 2)
        call my_equiv_3_2(xfd%VA, xfd%UNEW, (/ 1, 1, 1 /), (/ 1, 1 /),
     &                   (/ IV, 2 /), 2)
        call my_equiv_3_2(xfd%VA, xfd%U_AC, (/ 1, 1, IVX /), (/ 1, 1 /),
     &  
     &                   (/ IV, 1 /), 2)
        call my_equiv_3_2(xfd%VA, xfd%U_AC, (/ 1, 1, IVX /), (/ 1, 1 /),
     &  
     &                   (/ IV, 2 /), 2)
   83 CONTINUE
C
C     DP mod: copy from VB to QNEW to remove need for EQUIVALENCE
C       and from VB to Q_AC
      DO 84 IV = 1, IQX
        call my_equiv_3_1(xfd%VB, xfd%QNEW, (/ 1, 1, 1 /), 1, IV, 2)
        call my_equiv_3_1(xfd%VB, xfd%Q_AC, (/ 1, 1, IVX /), 1, IV, 2)
C
   84 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Adds on Newton deltas to boundary layer variables.
C      Checks for excessive changes and underrelaxes if necessary.
C      Calculates max and rms changes.
C      Also calculates the change in the global variable "AC".
C        If LALFA=.TRUE. , "AC" is CL
C        If LALFA=.FALSE., "AC" is alpha
C
C===================================================================70
      SUBROUTINE UPDATE(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C      REAL*8 :: UNEW(IVX,2), U_AC(IVX,2)
C      REAL*8 :: QNEW(IQX),   Q_AC(IQX)
C      EQUIVALENCE (VA(1,1,1), UNEW(1,1)) ,
C     &            (VB(1,1,1), QNEW(1)  )
C      EQUIVALENCE (VA(1,1,IVX), U_AC(1,1)) ,
C     &            (VB(1,1,IVX), Q_AC(1)  )
      REAL*8 MSQ
C
C---- max allowable alpha changes per iteration
      DALMAX =  0.5*xfd%DTOR
      DALMIN = -0.5*xfd%DTOR
C
C---- max allowable CL change per iteration
      DCLMAX =  0.5
      DCLMIN = -0.5
      IF(xfd%MATYP.NE.1) DCLMIN = MAX(-0.5 , -0.9*xfd%CL)
C
      HSTINV = xfd%GAMM1*(xfd%MINF/xfd%QINF)**2 / (1.0 + 0.5*xfd%GAMM1
     &  *xfd%MINF**2)
C
C---- calculate new Ue distribution assuming no under-relaxation
C-    also set the sensitivity of Ue wrt to alpha or Re
      DO 1 IS=1, 2
        DO 10 IBL=2, xfd%NBL(IS)
          I = xfd%IPAN(IBL,IS)
C
          DUI    = 0.
          DUI_AC = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, xfd%NBL(JS)
              J  = xfd%IPAN(JBL,JS)
              JV = xfd%ISYS(JBL,JS)
              UE_M = -xfd%VTI(IBL,IS)*xfd%VTI(JBL,JS)*xfd%DIJ(I,J)
              DUI    = DUI    + UE_M*(xfd%MASS(JBL,JS)+xfd%VDEL(3,1,JV))
     &  
              DUI_AC = DUI_AC + UE_M*(            -xfd%VDEL(3,2,JV))
 1000       CONTINUE
  100     CONTINUE
C
C-------- UINV depends on "AC" only if "AC" is alpha
          IF(xfd%LALFA) THEN
           UINV_AC = 0.
          ELSE
           UINV_AC = xfd%UINV_A(IBL,IS)
          ENDIF
C
          xfd%UNEW(IBL,IS) = xfd%UINV(IBL,IS) + DUI
          xfd%U_AC(IBL,IS) = UINV_AC      + DUI_AC
C
   10   CONTINUE
    1 CONTINUE
C
C---- set new Qtan from new Ue with appropriate sign change
      DO 2 IS=1, 2
        DO 20 IBL=2, xfd%IBLTE(IS)
          I = xfd%IPAN(IBL,IS)
          xfd%QNEW(I) = xfd%VTI(IBL,IS)*xfd%UNEW(IBL,IS)
          xfd%Q_AC(I) = xfd%VTI(IBL,IS)*xfd%U_AC(IBL,IS)
   20   CONTINUE
    2 CONTINUE
C
C---- calculate new CL from this new Qtan
      SA = SIN(xfd%ALFA)
      CA = COS(xfd%ALFA)
C
      BETA = SQRT(1.0 - xfd%MINF**2)
      BETA_MSQ = -0.5/BETA
C
      BFAC     = 0.5*xfd%MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA)
     &         - BFAC        / (1.0 + BETA) * BETA_MSQ
C
      CLNEW = 0.
      CL_A  = 0.
      CL_MS = 0.
      CL_AC = 0.
C
      I = 1
      CGINC = 1.0 - (xfd%QNEW(I)/xfd%QINF)**2
      CPG1  = CGINC / (BETA + BFAC*CGINC)
      CPG1_MS = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_Q = -2.0*xfd%QNEW(I)/xfd%QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_AC = CPC_CPI*CPI_Q*xfd%Q_AC(I)
C
      DO 3 I=1, xfd%N
        IP = I+1
        IF(I.EQ.xfd%N) IP = 1
C
        CGINC = 1.0 - (xfd%QNEW(IP)/xfd%QINF)**2
        CPG2  = CGINC / (BETA + BFAC*CGINC)
        CPG2_MS = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_Q = -2.0*xfd%QNEW(IP)/xfd%QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_AC = CPC_CPI*CPI_Q*xfd%Q_AC(IP)
C
        DX   =  (xfd%X(IP) - xfd%X(I))*CA + (xfd%Y(IP) - xfd%Y(I))*SA
        DX_A = -(xfd%X(IP) - xfd%X(I))*SA + (xfd%Y(IP) - xfd%Y(I))*CA
C
        AG    = 0.5*(CPG2    + CPG1   )
        AG_MS = 0.5*(CPG2_MS + CPG1_MS)
        AG_AC = 0.5*(CPG2_AC + CPG1_AC)
C
        CLNEW = CLNEW + DX  *AG
        CL_A  = CL_A  + DX_A*AG
        CL_MS = CL_MS + DX  *AG_MS
        CL_AC = CL_AC + DX  *AG_AC
C
        CPG1    = CPG2
        CPG1_MS = CPG2_MS
        CPG1_AC = CPG2_AC
    3 CONTINUE
C
C---- initialize under-relaxation factor
      xfd%RLX = 1.0
C
      IF(xfd%LALFA) THEN
C===== alpha is prescribed: AC is CL
C
C----- set change in Re to account for CL changing, since Re = Re(CL)
       DAC = (CLNEW - xfd%CL) / (1.0 - CL_AC - CL_MS*2.0*xfd%MINF
     &  *xfd%MINF_CL)
C
C----- set under-relaxation factor if Re change is too large
       IF(xfd%RLX*DAC .GT. DCLMAX) xfd%RLX = DCLMAX/DAC
       IF(xfd%RLX*DAC .LT. DCLMIN) xfd%RLX = DCLMIN/DAC
C
      ELSE
C===== CL is prescribed: AC is alpha
C
C----- set change in alpha to drive CL to prescribed value
       DAC = (CLNEW - xfd%CLSPEC) / (0.0 - CL_AC - CL_A)
C
C----- set under-relaxation factor if alpha change is too large
       IF(xfd%RLX*DAC .GT. DALMAX) xfd%RLX = DALMAX/DAC
       IF(xfd%RLX*DAC .LT. DALMIN) xfd%RLX = DALMIN/DAC
C
      ENDIF
C
      xfd%RMSBL = 0.
      xfd%RMXBL = 0.
C
      DHI = 1.5
      DLO = -.5
C
C---- calculate changes in BL variables and under-relaxation if needed
      DO 4 IS=1, 2
        DO 40 IBL=2, xfd%NBL(IS)
          IV = xfd%ISYS(IBL,IS)
C


C-------- set changes without underrelaxation
          DCTAU = xfd%VDEL(1,1,IV) - DAC*xfd%VDEL(1,2,IV)
          DTHET = xfd%VDEL(2,1,IV) - DAC*xfd%VDEL(2,2,IV)
          DMASS = xfd%VDEL(3,1,IV) - DAC*xfd%VDEL(3,2,IV)
          DUEDG = xfd%UNEW(IBL,IS) + DAC*xfd%U_AC(IBL,IS)  - 
     &   xfd%UEDG(IBL,IS)
          DDSTR = (DMASS - xfd%DSTR(IBL,IS)*DUEDG)/xfd%UEDG(IBL,IS)
C
C-------- normalize changes
          IF(IBL.LT.xfd%ITRAN(IS)) DN1 = DCTAU / 10.0
          IF(IBL.GE.xfd%ITRAN(IS)) DN1 = DCTAU / xfd%CTAU(IBL,IS)
          DN2 = DTHET / xfd%THET(IBL,IS)
          DN3 = DDSTR / xfd%DSTR(IBL,IS)
          DN4 = ABS(DUEDG)/0.25
C
C-------- accumulate for rms change
          xfd%RMSBL = xfd%RMSBL + DN1**2 + DN2**2 + DN3**2 + DN4**2
C          
C-------- see if Ctau needs underrelaxation
          RDN1 = xfd%RLX*DN1
          IF(ABS(DN1) .GT. ABS(xfd%RMXBL)) THEN
           xfd%RMXBL = DN1
           IF(IBL.LT.xfd%ITRAN(IS)) xfd%VMXBL = 'n'
           IF(IBL.GE.xfd%ITRAN(IS)) xfd%VMXBL = 'C'
           xfd%IMXBL = IBL
           xfd%ISMXBL = IS
          ENDIF
          IF(RDN1 .GT. DHI) xfd%RLX = DHI/DN1
          IF(RDN1 .LT. DLO) xfd%RLX = DLO/DN1
C
C-------- see if Theta needs underrelaxation
          RDN2 = xfd%RLX*DN2
          IF(ABS(DN2) .GT. ABS(xfd%RMXBL)) THEN
           xfd%RMXBL = DN2
           xfd%VMXBL = 'T'
           xfd%IMXBL = IBL
           xfd%ISMXBL = IS
          ENDIF
          IF(RDN2 .GT. DHI) xfd%RLX = DHI/DN2
          IF(RDN2 .LT. DLO) xfd%RLX = DLO/DN2
C
C-------- see if Dstar needs underrelaxation
          RDN3 = xfd%RLX*DN3
          IF(ABS(DN3) .GT. ABS(xfd%RMXBL)) THEN
           xfd%RMXBL = DN3
           xfd%VMXBL = 'D'
           xfd%IMXBL = IBL
           xfd%ISMXBL = IS
          ENDIF
          IF(RDN3 .GT. DHI) xfd%RLX = DHI/DN3
          IF(RDN3 .LT. DLO) xfd%RLX = DLO/DN3
C
C-------- see if Ue needs underrelaxation
          RDN4 = xfd%RLX*DN4
          IF(ABS(DN4) .GT. ABS(xfd%RMXBL)) THEN
           xfd%RMXBL = DUEDG
           xfd%VMXBL = 'U'
           xfd%IMXBL = IBL
           xfd%ISMXBL = IS
          ENDIF
          IF(RDN4 .GT. DHI) xfd%RLX = DHI/DN4
          IF(RDN4 .LT. DLO) xfd%RLX = DLO/DN4
C
   40   CONTINUE
    4 CONTINUE
C
C---- set true rms change
      xfd%RMSBL = SQRT( xfd%RMSBL / (4.0*FLOAT( xfd%NBL(1)+xfd%NBL(2) ))
     &   )
C
C
      IF(xfd%LALFA) THEN
C----- set underrelaxed change in Reynolds number from change in lift
       xfd%CL = xfd%CL + xfd%RLX*DAC
      ELSE
C----- set underrelaxed change in alpha
       xfd%ALFA = xfd%ALFA + xfd%RLX*DAC
       xfd%ADEG = xfd%ALFA/xfd%DTOR
      ENDIF
C
C---- update BL variables with underrelaxed changes
      DO 5 IS=1, 2
        DO 50 IBL=2, xfd%NBL(IS)
          IV = xfd%ISYS(IBL,IS)
C
          DCTAU = xfd%VDEL(1,1,IV) - DAC*xfd%VDEL(1,2,IV)
          DTHET = xfd%VDEL(2,1,IV) - DAC*xfd%VDEL(2,2,IV)
          DMASS = xfd%VDEL(3,1,IV) - DAC*xfd%VDEL(3,2,IV)
          DUEDG = xfd%UNEW(IBL,IS) + DAC*xfd%U_AC(IBL,IS)  - 
     &   xfd%UEDG(IBL,IS)
          DDSTR = (DMASS - xfd%DSTR(IBL,IS)*DUEDG)/xfd%UEDG(IBL,IS)
C
          xfd%CTAU(IBL,IS) = xfd%CTAU(IBL,IS) + xfd%RLX*DCTAU
          xfd%THET(IBL,IS) = xfd%THET(IBL,IS) + xfd%RLX*DTHET
          xfd%DSTR(IBL,IS) = xfd%DSTR(IBL,IS) + xfd%RLX*DDSTR
          xfd%UEDG(IBL,IS) = xfd%UEDG(IBL,IS) + xfd%RLX*DUEDG
C
          IF(IBL.GT.xfd%IBLTE(IS)) THEN
           IW = IBL - xfd%IBLTE(IS)
           DSWAKI = xfd%WGAP(IW)
          ELSE
           DSWAKI = 0.
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.xfd%ITRAN(IS))
     &      xfd%CTAU(IBL,IS) = MIN( xfd%CTAU(IBL,IS) , 0.25 )
C
          IF(IBL.LE.xfd%IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = xfd%UEDG(IBL,IS)**2*HSTINV
     &        / (xfd%GAMM1*(1.0 - 0.5*xfd%UEDG(IBL,IS)**2*HSTINV))
          DSW = xfd%DSTR(IBL,IS) - DSWAKI
          CALL DSLIM(DSW,xfd%THET(IBL,IS),xfd%UEDG(IBL,IS),MSQ,HKLIM)
          xfd%DSTR(IBL,IS) = DSW + DSWAKI
C
C-------- set new mass defect (nonlinear update)
          xfd%MASS(IBL,IS) = xfd%DSTR(IBL,IS) * xfd%UEDG(IBL,IS)
C
   50   CONTINUE
C
C------ make sure there are no "islands" of negative Ue
        DO IBL = 3, xfd%IBLTE(IS)
          IF(xfd%UEDG(IBL-1,IS) .GT. 0.0 .AND.
     &       xfd%UEDG(IBL  ,IS) .LE. 0.0       ) THEN
           xfd%UEDG(IBL,IS) = xfd%UEDG(IBL-1,IS)
           xfd%MASS(IBL,IS) = xfd%DSTR(IBL,IS) * xfd%UEDG(IBL,IS)
          ENDIF
        ENDDO
    5 CONTINUE
C
C
C---- equate upper wake arrays to lower wake arrays
      DO 6 KBL=1, xfd%NBL(2)-xfd%IBLTE(2)
        xfd%CTAU(xfd%IBLTE(1)+KBL,1) = xfd%CTAU(xfd%IBLTE(2)+KBL,2)
        xfd%THET(xfd%IBLTE(1)+KBL,1) = xfd%THET(xfd%IBLTE(2)+KBL,2)
        xfd%DSTR(xfd%IBLTE(1)+KBL,1) = xfd%DSTR(xfd%IBLTE(2)+KBL,2)
        xfd%UEDG(xfd%IBLTE(1)+KBL,1) = xfd%UEDG(xfd%IBLTE(2)+KBL,2)
         xfd%TAU(xfd%IBLTE(1)+KBL,1) =  xfd%TAU(xfd%IBLTE(2)+KBL,2)
         xfd%DIS(xfd%IBLTE(1)+KBL,1) =  xfd%DIS(xfd%IBLTE(2)+KBL,2)
         xfd%CTQ(xfd%IBLTE(1)+KBL,1) =  xfd%CTQ(xfd%IBLTE(2)+KBL,2)
        xfd%DELT(xfd%IBLTE(1)+KBL,1) = xfd%DELT(xfd%IBLTE(2)+KBL,2)
        xfd%TSTR(xfd%IBLTE(1)+KBL,1) = xfd%TSTR(xfd%IBLTE(2)+KBL,2)
    6 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Copies "2" variables to "1" variables to remove need for
C      EQUIVALENCE and COMMON blocks
C
C===================================================================70
      subroutine com2_to_com1(xbd)

      use xfoil_data_mod
      implicit none
      type(xbl_data_type), intent(inout) :: xbd

      xbd%X1     = xbd%X2    
      xbd%U1     = xbd%U2    
      xbd%T1     = xbd%T2    
      xbd%D1     = xbd%D2    
      xbd%S1     = xbd%S2    
      xbd%AMPL1  = xbd%AMPL2 
      xbd%U1_UEI = xbd%U2_UEI 
      xbd%U1_MS  = xbd%U2_MS  
      xbd%DW1    = xbd%DW2   
      xbd%H1     = xbd%H2    
      xbd%H1_T1  = xbd%H2_T2 
      xbd%H1_D1  = xbd%H2_D2 
      xbd%M1     = xbd%M2    
      xbd%M1_U1  = xbd%M2_U2 
      xbd%M1_MS  = xbd%M2_MS 
      xbd%R1     = xbd%R2    
      xbd%R1_U1  = xbd%R2_U2 
      xbd%R1_MS  = xbd%R2_MS 
      xbd%V1     = xbd%V2    
      xbd%V1_U1  = xbd%V2_U2  
      xbd%V1_MS  = xbd%V2_MS  
      xbd%V1_RE  = xbd%V2_RE  
      xbd%HK1    = xbd%HK2    
      xbd%HK1_U1 = xbd%HK2_U2  
      xbd%HK1_T1 = xbd%HK2_T2  
      xbd%HK1_D1 = xbd%HK2_D2  
      xbd%HK1_MS = xbd%HK2_MS  
      xbd%HS1    = xbd%HS2     
      xbd%HS1_U1 = xbd%HS2_U2    
      xbd%HS1_T1 = xbd%HS2_T2   
      xbd%HS1_D1 = xbd%HS2_D2  
      xbd%HS1_MS = xbd%HS2_MS  
      xbd%HS1_RE = xbd%HS2_RE  
      xbd%HC1    = xbd%HC2    
      xbd%HC1_U1 = xbd%HC2_U2   
      xbd%HC1_T1 = xbd%HC2_T2   
      xbd%HC1_D1 = xbd%HC2_D2   
      xbd%HC1_MS = xbd%HC2_MS   
      xbd%RT1    = xbd%RT2    
      xbd%RT1_U1 = xbd%RT2_U2  
      xbd%RT1_T1 = xbd%RT2_T2  
      xbd%RT1_MS = xbd%RT2_MS  
      xbd%RT1_RE = xbd%RT2_RE  
      xbd%CF1    = xbd%CF2    
      xbd%CF1_U1 = xbd%CF2_U2  
      xbd%CF1_T1 = xbd%CF2_T2  
      xbd%CF1_D1 = xbd%CF2_D2  
      xbd%CF1_MS = xbd%CF2_MS  
      xbd%CF1_RE = xbd%CF2_RE  
      xbd%DI1    = xbd%DI2    
      xbd%DI1_U1 = xbd%DI2_U2  
      xbd%DI1_T1 = xbd%DI2_T2  
      xbd%DI1_D1 = xbd%DI2_D2  
      xbd%DI1_S1 = xbd%DI2_S2  
      xbd%DI1_MS = xbd%DI2_MS  
      xbd%DI1_RE = xbd%DI2_RE  
      xbd%US1    = xbd%US2    
      xbd%US1_U1 = xbd%US2_U2  
      xbd%US1_T1 = xbd%US2_T2  
      xbd%US1_D1 = xbd%US2_D2  
      xbd%US1_MS = xbd%US2_MS  
      xbd%US1_RE = xbd%US2_RE  
      xbd%CQ1    = xbd%CQ2    
      xbd%CQ1_U1 = xbd%CQ2_U2  
      xbd%CQ1_T1 = xbd%CQ2_T2  
      xbd%CQ1_D1 = xbd%CQ2_D2  
      xbd%CQ1_MS = xbd%CQ2_MS  
      xbd%CQ1_RE = xbd%CQ2_RE  
      xbd%DE1    = xbd%DE2    
      xbd%DE1_U1 = xbd%DE2_U2  
      xbd%DE1_T1 = xbd%DE2_T2  
      xbd%DE1_D1 = xbd%DE2_D2  
      xbd%DE1_MS = xbd%DE2_MS  

      end subroutine com2_to_com1

C===================================================================70
C
C      Copies "1" variables to C1SAV array
C
C===================================================================70
      subroutine store_c1sav(xbd)

      use xfoil_data_mod
      implicit none
      type(xbl_data_type), intent(inout) :: xbd

      xbd%C1SAV(1)  = xbd%X1
      xbd%C1SAV(2)  = xbd%U1
      xbd%C1SAV(3)  = xbd%T1
      xbd%C1SAV(4)  = xbd%D1
      xbd%C1SAV(5)  = xbd%S1
      xbd%C1SAV(6)  = xbd%AMPL1
      xbd%C1SAV(7)  = xbd%U1_UEI
      xbd%C1SAV(8)  = xbd%U1_MS
      xbd%C1SAV(9)  = xbd%DW1
      xbd%C1SAV(10) = xbd%H1
      xbd%C1SAV(11) = xbd%H1_T1
      xbd%C1SAV(12) = xbd%H1_D1
      xbd%C1SAV(13) = xbd%M1
      xbd%C1SAV(14) = xbd%M1_U1
      xbd%C1SAV(15) = xbd%M1_MS
      xbd%C1SAV(16) = xbd%R1
      xbd%C1SAV(17) = xbd%R1_U1
      xbd%C1SAV(18) = xbd%R1_MS
      xbd%C1SAV(19) = xbd%V1
      xbd%C1SAV(20) = xbd%V1_U1
      xbd%C1SAV(21) = xbd%V1_MS
      xbd%C1SAV(22) = xbd%V1_RE
      xbd%C1SAV(23) = xbd%HK1
      xbd%C1SAV(24) = xbd%HK1_U1
      xbd%C1SAV(25) = xbd%HK1_T1
      xbd%C1SAV(26) = xbd%HK1_D1
      xbd%C1SAV(27) = xbd%HK1_MS
      xbd%C1SAV(28) = xbd%HS1
      xbd%C1SAV(29) = xbd%HS1_U1
      xbd%C1SAV(30) = xbd%HS1_T1
      xbd%C1SAV(31) = xbd%HS1_D1
      xbd%C1SAV(32) = xbd%HS1_MS
      xbd%C1SAV(33) = xbd%HS1_RE
      xbd%C1SAV(34) = xbd%HC1
      xbd%C1SAV(35) = xbd%HC1_U1
      xbd%C1SAV(36) = xbd%HC1_T1
      xbd%C1SAV(37) = xbd%HC1_D1
      xbd%C1SAV(38) = xbd%HC1_MS
      xbd%C1SAV(39) = xbd%RT1
      xbd%C1SAV(40) = xbd%RT1_U1
      xbd%C1SAV(41) = xbd%RT1_T1
      xbd%C1SAV(42) = xbd%RT1_MS
      xbd%C1SAV(43) = xbd%RT1_RE
      xbd%C1SAV(44) = xbd%CF1
      xbd%C1SAV(45) = xbd%CF1_U1
      xbd%C1SAV(46) = xbd%CF1_T1
      xbd%C1SAV(47) = xbd%CF1_D1
      xbd%C1SAV(48) = xbd%CF1_MS
      xbd%C1SAV(49) = xbd%CF1_RE
      xbd%C1SAV(50) = xbd%DI1
      xbd%C1SAV(51) = xbd%DI1_U1
      xbd%C1SAV(52) = xbd%DI1_T1
      xbd%C1SAV(53) = xbd%DI1_D1
      xbd%C1SAV(54) = xbd%DI1_S1
      xbd%C1SAV(55) = xbd%DI1_MS
      xbd%C1SAV(56) = xbd%DI1_RE
      xbd%C1SAV(57) = xbd%US1
      xbd%C1SAV(58) = xbd%US1_U1
      xbd%C1SAV(59) = xbd%US1_T1
      xbd%C1SAV(60) = xbd%US1_D1
      xbd%C1SAV(61) = xbd%US1_MS
      xbd%C1SAV(62) = xbd%US1_RE
      xbd%C1SAV(63) = xbd%CQ1
      xbd%C1SAV(64) = xbd%CQ1_U1
      xbd%C1SAV(65) = xbd%CQ1_T1
      xbd%C1SAV(66) = xbd%CQ1_D1
      xbd%C1SAV(67) = xbd%CQ1_MS
      xbd%C1SAV(68) = xbd%CQ1_RE
      xbd%C1SAV(69) = xbd%DE1
      xbd%C1SAV(70) = xbd%DE1_U1
      xbd%C1SAV(71) = xbd%DE1_T1
      xbd%C1SAV(72) = xbd%DE1_D1
      xbd%C1SAV(73) = xbd%DE1_MS

      end subroutine store_c1sav

C===================================================================70
C
C      Copies "1" variables from C1SAV array
C
C===================================================================70
      subroutine from_c1sav(xbd)

      use xfoil_data_mod
      implicit none
      type(xbl_data_type), intent(inout) :: xbd

      xbd%X1     = xbd%C1SAV(1)  
      xbd%U1     = xbd%C1SAV(2)    
      xbd%T1     = xbd%C1SAV(3)    
      xbd%D1     = xbd%C1SAV(4)      
      xbd%S1     = xbd%C1SAV(5)    
      xbd%AMPL1  = xbd%C1SAV(6)   
      xbd%U1_UEI = xbd%C1SAV(7)  
      xbd%U1_MS  = xbd%C1SAV(8)  
      xbd%DW1    = xbd%C1SAV(9)  
      xbd%H1     = xbd%C1SAV(10)   
      xbd%H1_T1  = xbd%C1SAV(11)  
      xbd%H1_D1  = xbd%C1SAV(12)  
      xbd%M1     = xbd%C1SAV(13)  
      xbd%M1_U1  = xbd%C1SAV(14)   
      xbd%M1_MS  = xbd%C1SAV(15) 
      xbd%R1     = xbd%C1SAV(16)  
      xbd%R1_U1  = xbd%C1SAV(17)
      xbd%R1_MS  = xbd%C1SAV(18)
      xbd%V1     = xbd%C1SAV(19)   
      xbd%V1_U1  = xbd%C1SAV(20) 
      xbd%V1_MS  = xbd%C1SAV(21)    
      xbd%V1_RE  = xbd%C1SAV(22)    
      xbd%HK1    = xbd%C1SAV(23)   
      xbd%HK1_U1 = xbd%C1SAV(24)    
      xbd%HK1_T1 = xbd%C1SAV(25)      
      xbd%HK1_D1 = xbd%C1SAV(26)        
      xbd%HK1_MS = xbd%C1SAV(27)       
      xbd%HS1    = xbd%C1SAV(28)    
      xbd%HS1_U1 = xbd%C1SAV(29)     
      xbd%HS1_T1 = xbd%C1SAV(30)     
      xbd%HS1_D1 = xbd%C1SAV(31)   
      xbd%HS1_MS = xbd%C1SAV(32)     
      xbd%HS1_RE = xbd%C1SAV(33)    
      xbd%HC1    = xbd%C1SAV(34)  
      xbd%HC1_U1 = xbd%C1SAV(35)   
      xbd%HC1_T1 = xbd%C1SAV(36)     
      xbd%HC1_D1 = xbd%C1SAV(37) 
      xbd%HC1_MS = xbd%C1SAV(38) 
      xbd%RT1    = xbd%C1SAV(39)  
      xbd%RT1_U1 = xbd%C1SAV(40) 
      xbd%RT1_T1 = xbd%C1SAV(41) 
      xbd%RT1_MS = xbd%C1SAV(42) 
      xbd%RT1_RE = xbd%C1SAV(43) 
      xbd%CF1    = xbd%C1SAV(44) 
      xbd%CF1_U1 = xbd%C1SAV(45) 
      xbd%CF1_T1 = xbd%C1SAV(46) 
      xbd%CF1_D1 = xbd%C1SAV(47) 
      xbd%CF1_MS = xbd%C1SAV(48) 
      xbd%CF1_RE = xbd%C1SAV(49) 
      xbd%DI1    = xbd%C1SAV(50)   
      xbd%DI1_U1 = xbd%C1SAV(51) 
      xbd%DI1_T1 = xbd%C1SAV(52) 
      xbd%DI1_D1 = xbd%C1SAV(53) 
      xbd%DI1_S1 = xbd%C1SAV(54) 
      xbd%DI1_MS = xbd%C1SAV(55) 
      xbd%DI1_RE = xbd%C1SAV(56) 
      xbd%US1    = xbd%C1SAV(57)   
      xbd%US1_U1 = xbd%C1SAV(58) 
      xbd%US1_T1 = xbd%C1SAV(59) 
      xbd%US1_D1 = xbd%C1SAV(60) 
      xbd%US1_MS = xbd%C1SAV(61) 
      xbd%US1_RE = xbd%C1SAV(62) 
      xbd%CQ1    = xbd%C1SAV(63)  
      xbd%CQ1_U1 = xbd%C1SAV(64) 
      xbd%CQ1_T1 = xbd%C1SAV(65) 
      xbd%CQ1_D1 = xbd%C1SAV(66) 
      xbd%CQ1_MS = xbd%C1SAV(67) 
      xbd%CQ1_RE = xbd%C1SAV(68) 
      xbd%DE1    = xbd%C1SAV(69)  
      xbd%DE1_U1 = xbd%C1SAV(70) 
      xbd%DE1_T1 = xbd%C1SAV(71) 
      xbd%DE1_D1 = xbd%C1SAV(72) 
      xbd%DE1_MS = xbd%C1SAV(73) 

      end subroutine from_c1sav

C===================================================================70
C
C      Copies "2" variables to C2SAV array
C
C===================================================================70
      subroutine store_c2sav(xbd)

      use xfoil_data_mod
      implicit none
      type(xbl_data_type), intent(inout) :: xbd

      xbd%C2SAV(1)  = xbd%X2
      xbd%C2SAV(2)  = xbd%U2
      xbd%C2SAV(3)  = xbd%T2
      xbd%C2SAV(4)  = xbd%D2
      xbd%C2SAV(5)  = xbd%S2
      xbd%C2SAV(6)  = xbd%AMPL2
      xbd%C2SAV(7)  = xbd%U2_UEI
      xbd%C2SAV(8)  = xbd%U2_MS
      xbd%C2SAV(9)  = xbd%DW2
      xbd%C2SAV(10) = xbd%H2
      xbd%C2SAV(11) = xbd%H2_T2
      xbd%C2SAV(12) = xbd%H2_D2
      xbd%C2SAV(13) = xbd%M2
      xbd%C2SAV(14) = xbd%M2_U2
      xbd%C2SAV(15) = xbd%M2_MS
      xbd%C2SAV(16) = xbd%R2
      xbd%C2SAV(17) = xbd%R2_U2
      xbd%C2SAV(18) = xbd%R2_MS
      xbd%C2SAV(19) = xbd%V2
      xbd%C2SAV(20) = xbd%V2_U2
      xbd%C2SAV(21) = xbd%V2_MS
      xbd%C2SAV(22) = xbd%V2_RE
      xbd%C2SAV(23) = xbd%HK2
      xbd%C2SAV(24) = xbd%HK2_U2
      xbd%C2SAV(25) = xbd%HK2_T2
      xbd%C2SAV(26) = xbd%HK2_D2
      xbd%C2SAV(27) = xbd%HK2_MS
      xbd%C2SAV(28) = xbd%HS2
      xbd%C2SAV(29) = xbd%HS2_U2
      xbd%C2SAV(30) = xbd%HS2_T2
      xbd%C2SAV(31) = xbd%HS2_D2
      xbd%C2SAV(32) = xbd%HS2_MS
      xbd%C2SAV(33) = xbd%HS2_RE
      xbd%C2SAV(34) = xbd%HC2
      xbd%C2SAV(35) = xbd%HC2_U2
      xbd%C2SAV(36) = xbd%HC2_T2
      xbd%C2SAV(37) = xbd%HC2_D2
      xbd%C2SAV(38) = xbd%HC2_MS
      xbd%C2SAV(39) = xbd%RT2
      xbd%C2SAV(40) = xbd%RT2_U2
      xbd%C2SAV(41) = xbd%RT2_T2
      xbd%C2SAV(42) = xbd%RT2_MS
      xbd%C2SAV(43) = xbd%RT2_RE
      xbd%C2SAV(44) = xbd%CF2
      xbd%C2SAV(45) = xbd%CF2_U2
      xbd%C2SAV(46) = xbd%CF2_T2
      xbd%C2SAV(47) = xbd%CF2_D2
      xbd%C2SAV(48) = xbd%CF2_MS
      xbd%C2SAV(49) = xbd%CF2_RE
      xbd%C2SAV(50) = xbd%DI2
      xbd%C2SAV(51) = xbd%DI2_U2
      xbd%C2SAV(52) = xbd%DI2_T2
      xbd%C2SAV(53) = xbd%DI2_D2
      xbd%C2SAV(54) = xbd%DI2_S2
      xbd%C2SAV(55) = xbd%DI2_MS
      xbd%C2SAV(56) = xbd%DI2_RE
      xbd%C2SAV(57) = xbd%US2
      xbd%C2SAV(58) = xbd%US2_U2
      xbd%C2SAV(59) = xbd%US2_T2
      xbd%C2SAV(60) = xbd%US2_D2
      xbd%C2SAV(61) = xbd%US2_MS
      xbd%C2SAV(62) = xbd%US2_RE
      xbd%C2SAV(63) = xbd%CQ2
      xbd%C2SAV(64) = xbd%CQ2_U2
      xbd%C2SAV(65) = xbd%CQ2_T2
      xbd%C2SAV(66) = xbd%CQ2_D2
      xbd%C2SAV(67) = xbd%CQ2_MS
      xbd%C2SAV(68) = xbd%CQ2_RE
      xbd%C2SAV(69) = xbd%DE2
      xbd%C2SAV(70) = xbd%DE2_U2
      xbd%C2SAV(71) = xbd%DE2_T2
      xbd%C2SAV(72) = xbd%DE2_D2
      xbd%C2SAV(73) = xbd%DE2_MS

      end subroutine store_c2sav

C===================================================================70
C
C      Copies "2" variables from C2SAV array
C
C===================================================================70
      subroutine from_c2sav(xbd)

      use xfoil_data_mod
      implicit none
      type(xbl_data_type), intent(inout) :: xbd

      xbd%X2     = xbd%C2SAV(1) 
      xbd%U2     = xbd%C2SAV(2) 
      xbd%T2     = xbd%C2SAV(3) 
      xbd%D2     = xbd%C2SAV(4) 
      xbd%S2     = xbd%C2SAV(5) 
      xbd%AMPL2  = xbd%C2SAV(6) 
      xbd%U2_UEI = xbd%C2SAV(7) 
      xbd%U2_MS  = xbd%C2SAV(8) 
      xbd%DW2    = xbd%C2SAV(9) 
      xbd%H2     = xbd%C2SAV(10)
      xbd%H2_T2  = xbd%C2SAV(11)
      xbd%H2_D2  = xbd%C2SAV(12)
      xbd%M2     = xbd%C2SAV(13)
      xbd%M2_U2  = xbd%C2SAV(14)
      xbd%M2_MS  = xbd%C2SAV(15)
      xbd%R2     = xbd%C2SAV(16)
      xbd%R2_U2  = xbd%C2SAV(17)
      xbd%R2_MS  = xbd%C2SAV(18)
      xbd%V2     = xbd%C2SAV(19)
      xbd%V2_U2  = xbd%C2SAV(20)
      xbd%V2_MS  = xbd%C2SAV(21)
      xbd%V2_RE  = xbd%C2SAV(22)
      xbd%HK2    = xbd%C2SAV(23)
      xbd%HK2_U2 = xbd%C2SAV(24)
      xbd%HK2_T2 = xbd%C2SAV(25)
      xbd%HK2_D2 = xbd%C2SAV(26)
      xbd%HK2_MS = xbd%C2SAV(27)
      xbd%HS2    = xbd%C2SAV(28)
      xbd%HS2_U2 = xbd%C2SAV(29)
      xbd%HS2_T2 = xbd%C2SAV(30)
      xbd%HS2_D2 = xbd%C2SAV(31)
      xbd%HS2_MS = xbd%C2SAV(32)
      xbd%HS2_RE = xbd%C2SAV(33)
      xbd%HC2    = xbd%C2SAV(34)
      xbd%HC2_U2 = xbd%C2SAV(35)
      xbd%HC2_T2 = xbd%C2SAV(36)
      xbd%HC2_D2 = xbd%C2SAV(37)
      xbd%HC2_MS = xbd%C2SAV(38)
      xbd%RT2    = xbd%C2SAV(39)
      xbd%RT2_U2 = xbd%C2SAV(40)
      xbd%RT2_T2 = xbd%C2SAV(41)
      xbd%RT2_MS = xbd%C2SAV(42)
      xbd%RT2_RE = xbd%C2SAV(43)
      xbd%CF2    = xbd%C2SAV(44)
      xbd%CF2_U2 = xbd%C2SAV(45)
      xbd%CF2_T2 = xbd%C2SAV(46)
      xbd%CF2_D2 = xbd%C2SAV(47)
      xbd%CF2_MS = xbd%C2SAV(48)
      xbd%CF2_RE = xbd%C2SAV(49)
      xbd%DI2    = xbd%C2SAV(50)
      xbd%DI2_U2 = xbd%C2SAV(51)
      xbd%DI2_T2 = xbd%C2SAV(52)
      xbd%DI2_D2 = xbd%C2SAV(53)
      xbd%DI2_S2 = xbd%C2SAV(54)
      xbd%DI2_MS = xbd%C2SAV(55)
      xbd%DI2_RE = xbd%C2SAV(56)
      xbd%US2    = xbd%C2SAV(57)
      xbd%US2_U2 = xbd%C2SAV(58)
      xbd%US2_T2 = xbd%C2SAV(59)
      xbd%US2_D2 = xbd%C2SAV(60)
      xbd%US2_MS = xbd%C2SAV(61)
      xbd%US2_RE = xbd%C2SAV(62)
      xbd%CQ2    = xbd%C2SAV(63)
      xbd%CQ2_U2 = xbd%C2SAV(64)
      xbd%CQ2_T2 = xbd%C2SAV(65)
      xbd%CQ2_D2 = xbd%C2SAV(66)
      xbd%CQ2_MS = xbd%C2SAV(67)
      xbd%CQ2_RE = xbd%C2SAV(68)
      xbd%DE2    = xbd%C2SAV(69)
      xbd%DE2_U2 = xbd%C2SAV(70)
      xbd%DE2_T2 = xbd%C2SAV(71)
      xbd%DE2_D2 = xbd%C2SAV(72)
      xbd%DE2_MS = xbd%C2SAV(73)

      end subroutine from_c2sav
