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
C===================================================================70
      SUBROUTINE BLPINI(bld)

      use xfoil_data_mod
      type(blpar_data_type), intent(inout) :: bld
C     
      bld%SCCON = 5.6
      bld%GACON = 6.70
      bld%GBCON = 0.75
      bld%GCCON = 18.0
      bld%DLCON =  0.9
C
      bld%CTRCON = 1.8
      bld%CTRCEX = 3.3
C
      bld%DUXCON = 1.0
C
      bld%CTCON = 0.5/(bld%GACON**2 * bld%GBCON)
C
      bld%CFFAC = 1.0
C
      RETURN
      END

C===================================================================70
C
C     Calculates current streamfunction Psi at panel node or wake node
C     I due to freestream and all bound vorticity Gam on the airfoil. 
C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
C
C     If SIGLIN=True, then Psi includes the effects of the viscous
C     source distribution Sig and the sensitivity vector dPsi/dSig
C     (DZDM) is calculated.
C
C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
C     is calculated, where n is the normal motion of the jth node.
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C
C===================================================================70
      SUBROUTINE PSILIN(xfd,I,XI,YI,NXI,NYI,PSI,PSI_NI,GEOLIN,SIGLIN)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd

      REAL*8 NXO, NYO, NXP, NYP, NXI, NYI
      LOGICAL GEOLIN,SIGLIN

C
C---- distance tolerance for determining if two points are the same
      SEPS = (xfd%S(xfd%N)-xfd%S(1)) * 1.0E-5
C
      IO = I
C
      xfd%COSA = COS(xfd%ALFA)
      xfd%SINA = SIN(xfd%ALFA)
C
      DO 3 JO=1, xfd%N
        xfd%DZDG(JO) = 0.0
        xfd%DZDN(JO) = 0.0
        xfd%DQDG(JO) = 0.0
    3 CONTINUE
C
      DO 4 JO=1, xfd%N
        xfd%DZDM(JO) = 0.0
        xfd%DQDM(JO) = 0.0
    4 CONTINUE
C
      xfd%Z_QINF = 0.
      xfd%Z_ALFA = 0.
      xfd%Z_QDOF0 = 0.
      xfd%Z_QDOF1 = 0.
      xfd%Z_QDOF2 = 0.
      xfd%Z_QDOF3 = 0.
C
      PSI    = 0.
      PSI_NI = 0.
C
      xfd%QTAN1 = 0.
      xfd%QTAN2 = 0.
      QTANM = 0.
C
      IF(xfd%SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = xfd%ANTE/xfd%DSTE
       SDS = xfd%ASTE/xfd%DSTE
      ENDIF
C
      DO 10 JO=1, xfd%N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.EQ.1) THEN
         JM = JO
        ELSE IF(JO.EQ.xfd%N-1) THEN
         JQ = JP
        ELSE IF(JO.EQ.xfd%N) THEN
         JP = 1
         IF((xfd%X(JO)-xfd%X(JP))**2 + (xfd%Y(JO)-xfd%Y(JP))**2 .LT.
     &   SEPS**2) GO TO 12
        ENDIF
C
        DSO = SQRT((xfd%X(JO)-xfd%X(JP))**2 + (xfd%Y(JO)-xfd%Y(JP))**2)
C
C------ skip null panel
        IF(DSO .EQ. 0.0) GO TO 10
C
        DSIO = 1.0 / DSO
C
        APAN = xfd%APANEL(JO)
C
        RX1 = XI - xfd%X(JO)
        RY1 = YI - xfd%Y(JO)
        RX2 = XI - xfd%X(JP)
        RY2 = YI - xfd%Y(JP)
C
        SX = (xfd%X(JP) - xfd%X(JO)) * DSIO
        SY = (xfd%Y(JP) - xfd%Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.xfd%N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        IF(IO.NE.JO .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*xfd%PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF(IO.NE.JP .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*xfd%PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = xfd%NX(JO)
         NYO = xfd%NY(JO)
         NXP = xfd%NX(JP)
         NYP = xfd%NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.EQ.xfd%N) GO TO 11
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*xfd%PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((xfd%X(JP)-xfd%X(JM))**2 + (xfd%Y(JP)-xfd%Y(JM))**2)
     &  
         DSIM = 1.0/DSM
C
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
C
         SSUM = (xfd%SIG(JP) - xfd%SIG(JO))*DSIO + (xfd%SIG(JP) -
     &   xfd%SIG(JM))*DSIM
         SDIF = (xfd%SIG(JP) - xfd%SIG(JO))*DSIO - (xfd%SIG(JP) -
     &   xfd%SIG(JM))*DSIM
C
         PSI = PSI + xfd%QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         xfd%DZDM(JM) = xfd%DZDM(JM) + xfd%QOPI*(-PSUM*DSIM + PDIF*DSIM)
     &  
         xfd%DZDM(JO) = xfd%DZDM(JO) + xfd%QOPI*(-PSUM*DSIO - PDIF*DSIO)
     &  
         xfd%DZDM(JP) = xfd%DZDM(JP) + xfd%QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + xfd%QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + xfd%QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         xfd%DQDM(JM) = xfd%DQDM(JM) + xfd%QOPI*(-PSNI*DSIM + PDNI*DSIM)
     &  
         xfd%DQDM(JO) = xfd%DQDM(JO) + xfd%QOPI*(-PSNI*DSIO - PDNI*DSIO)
     &  
         xfd%DQDM(JP) = xfd%DQDM(JP) + xfd%QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((xfd%X(JQ)-xfd%X(JO))**2 + (xfd%Y(JQ)-xfd%Y(JO))**2)
     &  
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (xfd%SIG(JQ) - xfd%SIG(JO))*DSIP + (xfd%SIG(JP) -
     &   xfd%SIG(JO))*DSIO
         SDIF = (xfd%SIG(JQ) - xfd%SIG(JO))*DSIP - (xfd%SIG(JP) -
     &   xfd%SIG(JO))*DSIO
C
         PSI = PSI + xfd%QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         xfd%DZDM(JO) = xfd%DZDM(JO) + xfd%QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         xfd%DZDM(JP) = xfd%DZDM(JP) + xfd%QOPI*( PSUM*DSIO - PDIF*DSIO)
     &  
         xfd%DZDM(JQ) = xfd%DZDM(JQ) + xfd%QOPI*( PSUM*DSIP + PDIF*DSIP)
     &  
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + xfd%QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         QTANM = QTANM + xfd%QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         xfd%DQDM(JO) = xfd%DQDM(JO) + xfd%QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         xfd%DQDM(JP) = xfd%DQDM(JP) + xfd%QOPI*( PSNI*DSIO - PDNI*DSIO)
     &  
         xfd%DQDM(JQ) = xfd%DQDM(JQ) + xfd%QOPI*( PSNI*DSIP + PDNI*DSIP)
     &  
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = xfd%GAMU(JP,1) + xfd%GAMU(JO,1)
        GSUM2 = xfd%GAMU(JP,2) + xfd%GAMU(JO,2)
        GDIF1 = xfd%GAMU(JP,1) - xfd%GAMU(JO,1)
        GDIF2 = xfd%GAMU(JP,2) - xfd%GAMU(JO,2)
C
        GSUM = xfd%GAM(JP) + xfd%GAM(JO)
        GDIF = xfd%GAM(JP) - xfd%GAM(JO)
C
        PSI = PSI + xfd%QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        xfd%DZDG(JO) = xfd%DZDG(JO) + xfd%QOPI*(PSIS-PSID)
        xfd%DZDG(JP) = xfd%DZDG(JP) + xfd%QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI + xfd%QOPI*(GSUM*PSNI + GDIF*PDNI)
C
        xfd%QTAN1 = xfd%QTAN1 + xfd%QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        xfd%QTAN2 = xfd%QTAN2 + xfd%QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        xfd%DQDG(JO) = xfd%DQDG(JO) + xfd%QOPI*(PSNI - PDNI)
        xfd%DQDG(JP) = xfd%DQDG(JP) + xfd%QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         xfd%DZDN(JO) = xfd%DZDN(JO)+ xfd%QOPI*GSUM*(PSX1*X1O + PSX2*X2O
     &   + PSYY*YYO)
     &                      + xfd%QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY
     &  *YYO)
         xfd%DZDN(JP) = xfd%DZDN(JP)+ xfd%QOPI*GSUM*(PSX1*X1P + PSX2*X2P
     &   + PSYY*YYP)
     &                      + xfd%QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY
     &  *YYP)
C------- dPsi/dP
         xfd%Z_QDOF0 = xfd%Z_QDOF0
     &           + xfd%QOPI*((PSIS-PSID)*xfd%QF0(JO) + (PSIS+PSID)
     &  *xfd%QF0(JP))
         xfd%Z_QDOF1 = xfd%Z_QDOF1
     &           + xfd%QOPI*((PSIS-PSID)*xfd%QF1(JO) + (PSIS+PSID)
     &  *xfd%QF1(JP))
         xfd%Z_QDOF2 = xfd%Z_QDOF2
     &           + xfd%QOPI*((PSIS-PSID)*xfd%QF2(JO) + (PSIS+PSID)
     &  *xfd%QF2(JP))
         xfd%Z_QDOF3 = xfd%Z_QDOF3
     &           + xfd%QOPI*((PSIS-PSID)*xfd%QF3(JO) + (PSIS+PSID)
     &  *xfd%QF3(JP))
        ENDIF
C
C
   10 CONTINUE
C
   11 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(xfd%GAMU(JP,1) - xfd%GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(xfd%GAMU(JP,2) - xfd%GAMU(JO,2))
      GAMTE1 = -.5*SDS*(xfd%GAMU(JP,1) - xfd%GAMU(JO,1))
      GAMTE2 = -.5*SDS*(xfd%GAMU(JP,2) - xfd%GAMU(JO,2))
C
      xfd%SIGTE = 0.5*SCS*(xfd%GAM(JP) - xfd%GAM(JO))
      xfd%GAMTE = -.5*SDS*(xfd%GAM(JP) - xfd%GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + xfd%HOPI*(PSIG*xfd%SIGTE + PGAM*xfd%GAMTE)
C
C---- dPsi/dGam
      xfd%DZDG(JO) = xfd%DZDG(JO) - xfd%HOPI*PSIG*SCS*0.5
      xfd%DZDG(JP) = xfd%DZDG(JP) + xfd%HOPI*PSIG*SCS*0.5
C
      xfd%DZDG(JO) = xfd%DZDG(JO) + xfd%HOPI*PGAM*SDS*0.5
      xfd%DZDG(JP) = xfd%DZDG(JP) - xfd%HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + xfd%HOPI*(PSIGNI*xfd%SIGTE + PGAMNI*xfd%GAMTE)
C
      xfd%QTAN1 = xfd%QTAN1 + xfd%HOPI*(PSIGNI*SIGTE1 + PGAMNI*GAMTE1)
      xfd%QTAN2 = xfd%QTAN2 + xfd%HOPI*(PSIGNI*SIGTE2 + PGAMNI*GAMTE2)
C
      xfd%DQDG(JO) = xfd%DQDG(JO) - xfd%HOPI*(PSIGNI*0.5*SCS - PGAMNI*0
     &  .5*SDS)
      xfd%DQDG(JP) = xfd%DQDG(JP) + xfd%HOPI*(PSIGNI*0.5*SCS - PGAMNI*0
     &  .5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       xfd%DZDN(JO) = xfd%DZDN(JO)
     &          + xfd%HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)
     &  *xfd%SIGTE
     &          + xfd%HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)
     &  *xfd%GAMTE
       xfd%DZDN(JP) = xfd%DZDN(JP)
     &          + xfd%HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)
     &  *xfd%SIGTE
     &          + xfd%HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)
     &  *xfd%GAMTE
C
C----- dPsi/dP
       xfd%Z_QDOF0 = xfd%Z_QDOF0 + xfd%HOPI*PSIG*0.5*(xfd%QF0(JP)
     &  -xfd%QF0(JO))*SCS
     &                   - xfd%HOPI*PGAM*0.5*(xfd%QF0(JP)-xfd%QF0(JO))
     &  *SDS
       xfd%Z_QDOF1 = xfd%Z_QDOF1 + xfd%HOPI*PSIG*0.5*(xfd%QF1(JP)
     &  -xfd%QF1(JO))*SCS
     &                   - xfd%HOPI*PGAM*0.5*(xfd%QF1(JP)-xfd%QF1(JO))
     &  *SDS
       xfd%Z_QDOF2 = xfd%Z_QDOF2 + xfd%HOPI*PSIG*0.5*(xfd%QF2(JP)
     &  -xfd%QF2(JO))*SCS
     &                   - xfd%HOPI*PGAM*0.5*(xfd%QF2(JP)-xfd%QF2(JO))
     &  *SDS
       xfd%Z_QDOF3 = xfd%Z_QDOF3 + xfd%HOPI*PSIG*0.5*(xfd%QF3(JP)
     &  -xfd%QF3(JO))*SCS
     &                   - xfd%HOPI*PGAM*0.5*(xfd%QF3(JP)-xfd%QF3(JO))
     &  *SDS
C
      ENDIF
C
   12 CONTINUE
C
C**** Freestream terms
      PSI = PSI + xfd%QINF*(xfd%COSA*YI - xfd%SINA*XI)
C
C---- dPsi/dn
      PSI_NI = PSI_NI + xfd%QINF*(xfd%COSA*NYI - xfd%SINA*NXI)
C
      xfd%QTAN1 = xfd%QTAN1 + xfd%QINF*NYI
      xfd%QTAN2 = xfd%QTAN2 - xfd%QINF*NXI
C
C---- dPsi/dQinf
      xfd%Z_QINF = xfd%Z_QINF + (xfd%COSA*YI - xfd%SINA*XI)
C
C---- dPsi/dalfa
      xfd%Z_ALFA = xfd%Z_ALFA - xfd%QINF*(xfd%SINA*YI + xfd%COSA*XI)
C
C     DP note: In Xfoil there is more stuff below here IF(LIMAGE), 
C     but LIMAGE is always .FALSE., even in Xfoil
C
      RETURN
      END !PSILIN

C===================================================================70
C
C     Calculates two surface vorticity (gamma) distributions
C     for alpha = 0, 90  degrees.  These are superimposed
C     in SPECAL or SPECCL for specified alpha or CL.
C
C===================================================================70
      SUBROUTINE GGCALC(xfd)

      use xfoil_data_mod
      use my_equivalence, only : my_equiv_3_2

      type(xfoil_data_type), intent(inout) :: xfd

C
C---- distance of internal control point ahead of sharp TE
C-    (fraction of smaller panel length adjacent to TE)
      BWT = 0.1
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE)
     &  WRITE(*,*) 'Calculating unit vorticity distributions ...'
C
      DO 10 I=1, xfd%N
        xfd%GAM(I) = 0.
        xfd%GAMU(I,1) = 0.
        xfd%GAMU(I,2) = 0.
   10 CONTINUE
      xfd%PSIO = 0.
C
C---- Set up matrix system for  Psi = Psio  on airfoil surface.
C-    The unknowns are (dGamma)i and dPsio.
      DO 20 I=1, xfd%N
C
C------ calculate Psi and dPsi/dGamma array for current node
        CALL PSILIN(xfd,I,xfd%X(I),xfd%Y(I),xfd%NX(I),xfd%NY(I),PSI,
     &  PSI_N,.FALSE.,.TRUE.)
C
        PSIINF = xfd%QINF*(COS(xfd%ALFA)*xfd%Y(I) - SIN(xfd%ALFA)
     &  *xfd%X(I))
C
C------ RES1 = PSI( 0) - PSIO
C------ RES2 = PSI(90) - PSIO
        RES1 =  xfd%QINF*xfd%Y(I)
        RES2 = -xfd%QINF*xfd%X(I)
C
C------ dRes/dGamma
        DO 201 J=1, xfd%N
          xfd%AIJ(I,J) = xfd%DZDG(J)
  201   CONTINUE
C
        DO 202 J=1, xfd%N
          xfd%BIJ(I,J) = -xfd%DZDM(J)

C         DP mod: copy to VM; used to replace equivalence statement
          call my_equiv_3_2(xfd%VM, xfd%BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     &    (/ I, J /), 1)

  202   CONTINUE
C
C------ dRes/dPsio
        xfd%AIJ(I,xfd%N+1) = -1.0
C
        xfd%GAMU(I,1) = -RES1
        xfd%GAMU(I,2) = -RES2
C
   20 CONTINUE
C
C---- set Kutta condition
C-    RES = GAM(1) + GAM(N)
      RES = 0.
C
      DO 30 J=1, xfd%N+1
        xfd%AIJ(xfd%N+1,J) = 0.0
   30 CONTINUE
C
      xfd%AIJ(xfd%N+1,1) = 1.0
      xfd%AIJ(xfd%N+1,xfd%N) = 1.0
C
      xfd%GAMU(xfd%N+1,1) = -RES
      xfd%GAMU(xfd%N+1,2) = -RES
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=1, xfd%N
        xfd%BIJ(xfd%N+1,J) = 0.

C       DP mod: copy to VM; used to replace equivalence statement
        call my_equiv_3_2(xfd%VM, xfd%BIJ, (/ 1, 1, 1 /), (/ 1, 1 /), 
     &                   (/ xfd%N+1, J /), 1)

   32 CONTINUE
C
      IF(xfd%SHARP) THEN
C----- set zero internal velocity in TE corner 
C
C----- set TE bisector angle
       AG1 = ATAN2(-xfd%YP(1),-xfd%XP(1)    )
       AG2 = ATANC( xfd%YP(xfd%N), xfd%XP(xfd%N),AG1)
       ABIS = 0.5*(AG1+AG2)
       CBIS = COS(ABIS)
       SBIS = SIN(ABIS)
C
C----- minimum panel length adjacent to TE
       DS1 = SQRT( (xfd%X(1)-xfd%X(2)  )**2 + (xfd%Y(1)-xfd%Y(2)  )**2 )
       DS2 = SQRT( (xfd%X(xfd%N)-xfd%X(xfd%N-1))**2 + (xfd%Y(xfd%N)
     &  - xfd%Y(xfd%N-1))**2 )
       DSMIN = MIN( DS1 , DS2 )
C
C----- control point on bisector just ahead of TE point
       XBIS = xfd%XTE - BWT*DSMIN*CBIS
       YBIS = xfd%YTE - BWT*DSMIN*SBIS
ccc       write(*,*) xbis, ybis
C
C----- set velocity component along bisector line
       CALL PSILIN(xfd,0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,.FALSE.,.TRUE.)
C
CCC--- RES = DQDGj*Gammaj + DQDMj*Massj + QINF*(COSA*CBIS + SINA*SBIS)
       RES = QBIS
C
C----- dRes/dGamma
       DO J=1, xfd%N
         xfd%AIJ(xfd%N,J) = xfd%DQDG(J)
       ENDDO
C
C----- -dRes/dMass
       DO J=1, xfd%N
         xfd%BIJ(xfd%N,J) = -xfd%DQDM(J)

C        DP mod: Copy to VM used to replace equivalence statement
         call my_equiv_3_2(xfd%VM, xfd%BIJ, (/ 1, 1, 1 /), (/ 1, 1 /), 
     &                    (/ xfd%N, J /), 1)

       ENDDO
C
C----- dRes/dPsio
       xfd%AIJ(xfd%N,xfd%N+1) = 0.
C
C----- -dRes/dUinf
       xfd%GAMU(xfd%N,1) = -CBIS
C
C----- -dRes/dVinf
       xfd%GAMU(xfd%N,2) = -SBIS
C
      ENDIF
C
C---- LU-factor coefficient matrix AIJ
      CALL LUDCMP(IQX,xfd%N+1,xfd%AIJ,xfd%AIJPIV)
      xfd%LQAIJ = .TRUE.
C
C---- solve system for the two vorticity distributions
      CALL BAKSUB(IQX,xfd%N+1,xfd%AIJ,xfd%AIJPIV,xfd%GAMU(1,1))
      CALL BAKSUB(IQX,xfd%N+1,xfd%AIJ,xfd%AIJPIV,xfd%GAMU(1,2))
C
C---- set inviscid alpha=0,90 surface speeds for this geometry
      DO 50 I=1, xfd%N
        xfd%QINVU(I,1) = xfd%GAMU(I,1)
        xfd%QINVU(I,2) = xfd%GAMU(I,2)
   50 CONTINUE
C
      xfd%LGAMU = .TRUE.

      RETURN
      END

C===================================================================70
C
C     Sets inviscid panel tangential velocity for
C     current alpha.
C
C===================================================================70
      SUBROUTINE QISET(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      xfd%COSA = COS(xfd%ALFA)
      xfd%SINA = SIN(xfd%ALFA)
C
      DO 5 I=1, xfd%N+xfd%NW
        xfd%QINV  (I) =  xfd%COSA*xfd%QINVU(I,1) + xfd%SINA*xfd%QINVU(I
     &  ,2)
        xfd%QINV_A(I) = -xfd%SINA*xfd%QINVU(I,1) + xfd%COSA*xfd%QINVU(I
     &  ,2)
    5 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Integrates surface pressures to get CL and CM.
C     Integrates skin friction to get CDF.
C     Calculates dCL/dAlpha for prescribed-CL routines.
C
C===================================================================70
      SUBROUTINE CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, 
     &                  XREF,YREF,
     &                  CL,CM,CDP, CL_ALF,CL_MSQ)
      DIMENSION X(N),Y(N), GAM(N), GAM_A(N)
      REAL*8 MINF
C
ccC---- moment-reference coordinates
cc      XREF = 0.25
cc      YREF = 0.
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      BFAC     = 0.5*MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA)
     &         - BFAC        / (1.0 + BETA) * BETA_MSQ
C
      CL = 0.0
      CM = 0.0

      CDP = 0.0
C
      CL_ALF = 0.
      CL_MSQ = 0.
C
      I = 1
      CGINC = 1.0 - (GAM(I)/QINF)**2
      CPG1     = CGINC/(BETA + BFAC*CGINC)
      CPG1_MSQ = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_GAM = -2.0*GAM(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_ALF = CPC_CPI*CPI_GAM*GAM_A(I)
C
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC = 1.0 - (GAM(IP)/QINF)**2
        CPG2     = CGINC/(BETA + BFAC*CGINC)
        CPG2_MSQ = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_GAM = -2.0*GAM(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_ALF = CPC_CPI*CPI_GAM*GAM_A(IP)
C
        DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DY = (Y(IP) - Y(I))*CA - (X(IP) - X(I))*SA
        DG = CPG2 - CPG1
C
        AX = (0.5*(X(IP)+X(I))-XREF)*CA + (0.5*(Y(IP)+Y(I))-YREF)*SA
        AY = (0.5*(Y(IP)+Y(I))-YREF)*CA - (0.5*(X(IP)+X(I))-XREF)*SA
        AG = 0.5*(CPG2 + CPG1)
C
        DX_ALF = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
        AG_ALF = 0.5*(CPG2_ALF + CPG1_ALF)
        AG_MSQ = 0.5*(CPG2_MSQ + CPG1_MSQ)
C
        CL     = CL     + DX* AG
        CDP    = CDP    - DY* AG
        CM     = CM     - DX*(AG*AX + DG*DX/12.0)
     &                  - DY*(AG*AY + DG*DY/12.0)
C
        CL_ALF = CL_ALF + DX*AG_ALF + AG*DX_ALF
        CL_MSQ = CL_MSQ + DX*AG_MSQ
C
        CPG1 = CPG2
        CPG1_ALF = CPG2_ALF
        CPG1_MSQ = CPG2_MSQ
   10 CONTINUE
C
      RETURN
      END ! CLCALC

C===================================================================70
C
C     Sets compressible Cp from speed.
C
C===================================================================70
      SUBROUTINE CPCALC(N,Q,QINF,MINF,CP,SILENT_MODE)

      use iso_c_binding

      DIMENSION Q(N),CP(N)
      REAL*8 MINF
      LOGICAL(c_bool) SILENT_MODE
C
      LOGICAL DENNEG
C
      BETA = SQRT(1.0 - MINF**2)
      BFAC = 0.5*MINF**2 / (1.0 + BETA)
C
      DENNEG = .FALSE.
C
      DO 20 I=1, N
        CPINC = 1.0 - (Q(I)/QINF)**2
        DEN = BETA + BFAC*CPINC
        CP(I) = CPINC / DEN
        IF(DEN .LE. 0.0) DENNEG = .TRUE.
  20  CONTINUE
C
C     DP mod: added SILENT_MODE option
      IF(DENNEG .AND. .NOT. SILENT_MODE) THEN
       WRITE(*,*)
       WRITE(*,*) 'CPCALC: Local speed too large. ',
     &            'Compressibility corrections invalid.'
      ENDIF
C
      RETURN
      END ! CPCALC

C===================================================================70
C
C Sets MINF from input
C
C===================================================================70
      SUBROUTINE MINFSET(xfd,MACH_INPUT)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd

      REAL*8 MACH_INPUT

      IF(xfd%MINF1.GE.1.0) THEN
        WRITE(*,*) 'Supersonic freestream not allowed'
        STOP
      ENDIF
      xfd%MINF1 = MACH_INPUT

      CALL MRCL(xfd,1.0,xfd%MINF_CL,REINF_CL)
      CALL COMSET(xfd)

      IF ( (xfd%MINF.GT.0.0) .AND. (.NOT. xfd%SILENT_MODE) ) 
     &   WRITE(*,1300) xfd%CPSTAR, xfd%QSTAR/xfd%QINF
 1300 FORMAT(/' Sonic Cp =', F10.2, '      Sonic Q/Qinf =', F10.3/)

      CALL CPCALC(xfd%N,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      IF(xfd%LVISC) CALL CPCALC(xfd%N+xfd%NW,xfd%QVIS,xfd%QINF,xfd%MINF
     &  ,xfd%CPV,xfd%SILENT_MODE)
      xfd%LVCONV = .FALSE.

      END ! MINFSET

C===================================================================70
C
C     Converges to specified alpha.
C
C===================================================================70
      SUBROUTINE SPECAL(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd

      REAL*8 MINF_CLM, MSQ_CLM

C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.xfd%LGAMU .OR. .NOT.xfd%LQAIJ) CALL GGCALC(xfd)
C
      xfd%COSA = COS(xfd%ALFA)
      xfd%SINA = SIN(xfd%ALFA)
C
C---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, xfd%N
        xfd%GAM(I)   =  xfd%COSA*xfd%GAMU(I,1) + xfd%SINA*xfd%GAMU(I,2)
        xfd%GAM_A(I) = -xfd%SINA*xfd%GAMU(I,1) + xfd%COSA*xfd%GAMU(I,2)
   50 CONTINUE
      xfd%PSIO = xfd%COSA*xfd%GAMU(xfd%N+1,1) + xfd%SINA*xfd%GAMU(xfd%N
     &  +1,2)
C
      CALL TECALC
      CALL QISET(xfd)
C
C---- set initial guess for the Newton variable CLM
      CLM = 1.0
C
C---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(xfd,CLM,MINF_CLM,REINF_CLM)
      CALL COMSET(xfd)
C
C---- set corresponding CL(M)
      CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA,xfd%MINF
     &  ,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &            xfd%CL,xfd%CM,xfd%CDP, xfd%CL_ALF,xfd%CL_MSQ)
C
C---- iterate on CLM
      DO 100 ITCL=1, 20
C
        MSQ_CLM = 2.0*xfd%MINF*MINF_CLM
        DCLM = (xfd%CL - CLM)/(1.0 - xfd%CL_MSQ*MSQ_CLM)
C
        CLM1 = CLM
        xfd%RLX = 1.0
C
C------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
C
          CLM = CLM1 + xfd%RLX*DCLM
C
C-------- set new freestream Mach M(CLM)
          CALL MRCL(xfd,CLM,MINF_CLM,REINF_CLM)
C
C-------- if Mach is OK, go do next Newton iteration
          IF(xfd%MATYP.EQ.1 .OR. xfd%MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0)
     &   GO TO 91
C
          xfd%RLX = 0.5*xfd%RLX
   90   CONTINUE
   91   CONTINUE
C
C------ set new CL(M)
        CALL COMSET(xfd)
        CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA
     &  ,xfd%MINF,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &              xfd%CL,xfd%CM,xfd%CDP,xfd%CL_ALF,xfd%CL_MSQ)
C
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
C
  100 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE) 
     &  WRITE(*,*) 'SPECAL:  Minf convergence failed'
  110 CONTINUE
C
C---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(xfd,xfd%CL,xfd%MINF_CL,REINF_CL)
      CALL COMSET(xfd)
      CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA,xfd%MINF
     &  ,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &            xfd%CL,xfd%CM,xfd%CDP, xfd%CL_ALF,xfd%CL_MSQ)
      CALL CPCALC(xfd%N,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      IF(xfd%LVISC) THEN
       CALL CPCALC(xfd%N+xfd%NW,xfd%QVIS,xfd%QINF,xfd%MINF,xfd%CPV
     &  ,xfd%SILENT_MODE)
       CALL CPCALC(xfd%N+xfd%NW,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      ELSE
       CALL CPCALC(xfd%N,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      ENDIF
C      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECAL

C===================================================================70
C
C     Converges to specified inviscid CL.
C
C===================================================================70
      SUBROUTINE SPECCL(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.xfd%LGAMU .OR. .NOT.xfd%LQAIJ) CALL GGCALC(xfd)
C
C---- set freestream Mach from specified CL -- Mach will be held fixed
      CALL MRCL(xfd,xfd%CLSPEC,xfd%MINF_CL,REINF_CL)
      CALL COMSET(xfd)
C
C---- current alpha is the initial guess for Newton variable ALFA
      xfd%COSA = COS(xfd%ALFA)
      xfd%SINA = SIN(xfd%ALFA)
      DO 10 I=1, xfd%N
        xfd%GAM(I)   =  xfd%COSA*xfd%GAMU(I,1) + xfd%SINA*xfd%GAMU(I,2)
        xfd%GAM_A(I) = -xfd%SINA*xfd%GAMU(I,1) + xfd%COSA*xfd%GAMU(I,2)
   10 CONTINUE
      xfd%PSIO = xfd%COSA*xfd%GAMU(xfd%N+1,1) + xfd%SINA*xfd%GAMU(xfd%N
     &  +1,2)
C
C---- get corresponding CL, CL_alpha, CL_Mach
      CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA,xfd%MINF
     &  ,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &            xfd%CL,xfd%CM,xfd%CDP, xfd%CL_ALF,xfd%CL_MSQ)
C
C---- Newton loop for alpha to get specified inviscid CL
      DO 100 ITAL=1, 20
C
        DALFA = (xfd%CLSPEC - xfd%CL) / xfd%CL_ALF
        xfd%RLX = 1.0
C
        xfd%ALFA = xfd%ALFA + xfd%RLX*DALFA
C
C------ set new surface speed distribution
        xfd%COSA = COS(xfd%ALFA)
        xfd%SINA = SIN(xfd%ALFA)
        DO 40 I=1, xfd%N
          xfd%GAM(I)   =  xfd%COSA*xfd%GAMU(I,1) + xfd%SINA*xfd%GAMU(I,2
     &  )
          xfd%GAM_A(I) = -xfd%SINA*xfd%GAMU(I,1) + xfd%COSA*xfd%GAMU(I,2
     &  )
   40   CONTINUE
        xfd%PSIO = xfd%COSA*xfd%GAMU(xfd%N+1,1) + xfd%SINA
     &  *xfd%GAMU(xfd%N+1,2)
C
C------ set new CL(alpha)
        CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA
     &  ,xfd%MINF,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &              xfd%CL,xfd%CM,xfd%CDP,xfd%CL_ALF,xfd%CL_MSQ)
C
        IF(ABS(DALFA).LE.1.0E-6) GO TO 110
  100 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE) 
     &  WRITE(*,*) 'SPECCL:  xfd%CL convergence failed'
  110 CONTINUE
C
C---- set final surface speed and Cp distributions
      CALL TECALC
      CALL QISET(xfd)
      IF(xfd%LVISC) THEN
       CALL CPCALC(xfd%N+xfd%NW,xfd%QVIS,xfd%QINF,xfd%MINF,xfd%CPV
     &  ,xfd%SILENT_MODE)
       CALL CPCALC(xfd%N+xfd%NW,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      ELSE
       CALL CPCALC(xfd%N,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      ENDIF
C      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECCL

C===================================================================70
C
C     Sets inviscid tangential velocity for alpha = 0, 90
C     on wake due to freestream and airfoil surface vorticity.
C
C===================================================================70
      SUBROUTINE QWCALC(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
C---- first wake point (same as TE)
      xfd%QINVU(xfd%N+1,1) = xfd%QINVU(xfd%N,1)
      xfd%QINVU(xfd%N+1,2) = xfd%QINVU(xfd%N,2)
C
C---- rest of wake
      DO 10 I=xfd%N+2, xfd%N+xfd%NW
        CALL PSILIN(xfd,I,xfd%X(I),xfd%Y(I),xfd%NX(I),xfd%NY(I),PSI,
     &  PSI_NI,.FALSE.,.FALSE.)
        xfd%QINVU(I,1) = xfd%QTAN1
        xfd%QINVU(I,2) = xfd%QTAN2
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE GAMQV(xfd)
      
      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      DO 10 I=1, xfd%N
        xfd%GAM(I)   = xfd%QVIS(I)
        xfd%GAM_A(I) = xfd%QINV_A(I)
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Locates stagnation point arc length 
C     location SST and panel index IST.
C
C===================================================================70
      SUBROUTINE STFIND(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      DO 10 I=1, xfd%N-1
        IF(xfd%GAM(I).GE.0.0 .AND. xfd%GAM(I+1).LT.0.0) GO TO 11
   10 CONTINUE
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE)
     &  WRITE(*,*) 'STFIND: Stagnation point not found. Continuing ...'
      I = xfd%N/2
C
   11 CONTINUE
C
      xfd%IST = I
      DGAM = xfd%GAM(I+1) - xfd%GAM(I)
      DS = xfd%S(I+1) - xfd%S(I)
C
C---- evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
      IF(xfd%GAM(I) .LT. -xfd%GAM(I+1)) THEN
       xfd%SST = xfd%S(I)   - DS*(xfd%GAM(I)  /DGAM)
      ELSE
       xfd%SST = xfd%S(I+1) - DS*(xfd%GAM(I+1)/DGAM)
      ENDIF
C
C---- tweak stagnation point if it falls right on a node (very unlikely)
      IF(xfd%SST .LE. xfd%S(I)  ) xfd%SST = xfd%S(I)   + 1.0E-7
      IF(xfd%SST .GE. xfd%S(I+1)) xfd%SST = xfd%S(I+1) - 1.0E-7
C
      xfd%SST_GO = (xfd%SST  - xfd%S(I+1))/DGAM
      xfd%SST_GP = (xfd%S(I) - xfd%SST   )/DGAM
C
      RETURN
      END

C===================================================================70
C
C     Sets  BL location -> panel location  pointer array IPAN
C
C===================================================================70
      SUBROUTINE IBLPAN(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
C---- top surface first
      IS = 1
C
      IBL = 1
      DO 10 I=xfd%IST, 1, -1
        IBL = IBL+1
        xfd%IPAN(IBL,IS) = I
        xfd%VTI(IBL,IS) = 1.0
   10 CONTINUE
C
      xfd%IBLTE(IS) = IBL
      xfd%NBL(IS) = IBL
C
C---- bottom surface next
      IS = 2
C
      IBL = 1
      DO 20 I=xfd%IST+1, xfd%N
        IBL = IBL+1
        xfd%IPAN(IBL,IS) = I
        xfd%VTI(IBL,IS) = -1.0
   20 CONTINUE
C
C---- wake
      xfd%IBLTE(IS) = IBL
C
      DO 25 IW=1, xfd%NW
        I = xfd%N+IW
        IBL = xfd%IBLTE(IS)+IW
        xfd%IPAN(IBL,IS) = I
         xfd%VTI(IBL,IS) = -1.0
   25 CONTINUE
C
      xfd%NBL(IS) = xfd%IBLTE(IS) + xfd%NW
C     DP note: stopped here
C
C---- upper wake pointers (for plotting only)
      DO 35 IW=1, xfd%NW
        xfd%IPAN(xfd%IBLTE(1)+IW,1) = xfd%IPAN(xfd%IBLTE(2)+IW,2)
         xfd%VTI(xfd%IBLTE(1)+IW,1) = 1.0
   35 CONTINUE
C
C
      IBLMAX = MAX(xfd%IBLTE(1),xfd%IBLTE(2)) + xfd%NW
      IF(IBLMAX.GT.IVX) THEN
        WRITE(*,*) ' ***  BL array overflow.'
        WRITE(*,*) ' ***  Increase IVX to at least', IBLMAX
        STOP
      ENDIF
C
      xfd%LIPAN = .TRUE.
      RETURN
      END

C===================================================================70
C
C     Sets BL arc length array on each airfoil side and wake
C
C===================================================================70
      SUBROUTINE XICALC(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd

      DATA XFEPS / 1.0E-7 /
C
C---- minimum xi node arc length near stagnation point
      XEPS = XFEPS*(xfd%S(xfd%N)-xfd%S(1))
C
      IS = 1
C
      xfd%XSSI(1,IS) = 0.
C
      DO 10 IBL=2, xfd%IBLTE(IS)
        I = xfd%IPAN(IBL,IS)
        xfd%XSSI(IBL,IS) = MAX( xfd%SST - xfd%S(I) , XEPS )
   10 CONTINUE
C
C
      IS = 2
C
      xfd%XSSI(1,IS) = 0.
C
      DO 20 IBL=2, xfd%IBLTE(IS)
        I = xfd%IPAN(IBL,IS)
        xfd%XSSI(IBL,IS) = MAX( xfd%S(I) - xfd%SST , XEPS )
   20 CONTINUE
C
C
      IS1 = 1
      IS2 = 2
C
      IBL1 = xfd%IBLTE(IS1) + 1
      xfd%XSSI(IBL1,IS1) = xfd%XSSI(IBL1-1,IS1)
C
      IBL2 = xfd%IBLTE(IS2) + 1
      xfd%XSSI(IBL2,IS2) = xfd%XSSI(IBL2-1,IS2)
C
      DO 25 IBL=xfd%IBLTE(IS)+2, xfd%NBL(IS)
        I = xfd%IPAN(IBL,IS)
        DXSSI = SQRT((xfd%X(I)-xfd%X(I-1))**2 + (xfd%Y(I)-xfd%Y(I-1))**2
     &  )
C
        IBL1 = xfd%IBLTE(IS1) + IBL - xfd%IBLTE(IS)
        IBL2 = xfd%IBLTE(IS2) + IBL - xfd%IBLTE(IS)
        xfd%XSSI(IBL1,IS1) = xfd%XSSI(IBL1-1,IS1) + DXSSI
        xfd%XSSI(IBL2,IS2) = xfd%XSSI(IBL2-1,IS2) + DXSSI
   25 CONTINUE
C
C---- trailing edge flap length to TE gap ratio
      TELRAT = 2.50
C
C---- set up parameters for TE flap cubics
C
ccc   DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)    !!! BUG  2/2/95
C
      CROSP = (xfd%XP(1)*xfd%YP(xfd%N) - xfd%YP(1)*xfd%XP(xfd%N))
     &      / SQRT(  (xfd%XP(1)**2 + xfd%YP(1)**2)
     &              *(xfd%XP(xfd%N)**2 + xfd%YP(xfd%N)**2) )
      DWDXTE = CROSP / SQRT(1.0 - CROSP**2)
C
C---- limit cubic to avoid absurd TE gap widths
      DWDXTE = MAX(DWDXTE,-3.0/TELRAT)
      DWDXTE = MIN(DWDXTE, 3.0/TELRAT)
C
      AA =  3.0 + TELRAT*DWDXTE
      BB = -2.0 - TELRAT*DWDXTE
C
      IF(xfd%SHARP) THEN
       DO 30 IW=1, xfd%NW
         xfd%WGAP(IW) = 0.
   30  CONTINUE
      ELSE
C----- set TE flap (wake gap) array
       IS = 2
       DO 35 IW=1, xfd%NW
         IBL = xfd%IBLTE(IS) + IW
         ZN = 1.0 - (xfd%XSSI(IBL,IS)-xfd%XSSI(xfd%IBLTE(IS),IS)) /
     &   (TELRAT*xfd%ANTE)
         xfd%WGAP(IW) = 0.
         IF(ZN.GE.0.0) xfd%WGAP(IW) = xfd%ANTE * (AA + BB*ZN)*ZN**2
   35  CONTINUE
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Sets inviscid Ue from panel inviscid tangential velocity
C
C===================================================================70
      SUBROUTINE UICALC(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      DO 10 IS=1, 2
        xfd%UINV  (1,IS) = 0.
        xfd%UINV_A(1,IS) = 0.
        DO 110 IBL=2, xfd%NBL(IS)
          I = xfd%IPAN(IBL,IS)
          xfd%UINV  (IBL,IS) = xfd%VTI(IBL,IS)*xfd%QINV  (I)
          xfd%UINV_A(IBL,IS) = xfd%VTI(IBL,IS)*xfd%QINV_A(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets panel viscous tangential velocity from viscous Ue
C
C===================================================================70
      SUBROUTINE QVFUE(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      DO 1 IS=1, 2
        DO 10 IBL=2, xfd%NBL(IS)
          I = xfd%IPAN(IBL,IS)
          xfd%QVIS(I) = xfd%VTI(IBL,IS)*xfd%UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE CDCALC(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
      SA = SIN(xfd%ALFA)
      CA = COS(xfd%ALFA)
C
      IF(xfd%LVISC .AND. xfd%LBLINI) THEN
C
C----- set variables at the end of the wake
       THWAKE = xfd%THET(xfd%NBL(2),2)
       URAT   = xfd%UEDG(xfd%NBL(2),2)/xfd%QINF
       UEWAKE = xfd%UEDG(xfd%NBL(2),2) * (1.0-xfd%TKLAM) / (1.0 -
     &   xfd%TKLAM*URAT**2)
       SHWAKE = xfd%DSTR(xfd%NBL(2),2)/xfd%THET(xfd%NBL(2),2)
C
C----- extrapolate wake to downstream infinity using Squire-Young relation
C      (reduces errors of the wake not being long enough)
       xfd%CD = 2.0*THWAKE * (UEWAKE/xfd%QINF)**(0.5*(5.0+SHWAKE))
C
      ELSE
C
       xfd%CD = 0.0
C
      ENDIF
C
C---- calculate friction drag coefficient
      xfd%CDF = 0.0
      DO 20 IS=1, 2
        DO 205 IBL=3, xfd%IBLTE(IS)
          I  = xfd%IPAN(IBL  ,IS)
          IM = xfd%IPAN(IBL-1,IS)
          DX = (xfd%X(I) - xfd%X(IM))*CA + (xfd%Y(I) - xfd%Y(IM))*SA
          xfd%CDF = xfd%CDF + 0.5*(xfd%TAU(IBL,IS)+xfd%TAU(IBL-1,IS))*DX
     &   * 2.0/xfd%QINF**2
 205    CONTINUE
 20   CONTINUE
C
      RETURN
      END ! CDCALC

C===================================================================70
C
C     Moves stagnation point location to new panel.
C
C===================================================================70
      SUBROUTINE STMOVE(xfd)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
C
C---- locate new stagnation point arc length SST from GAM distribution
      ISTOLD = xfd%IST
      CALL STFIND(xfd)
C
      IF(ISTOLD.EQ.xfd%IST) THEN
C
C----- recalculate new arc length array
       CALL XICALC(xfd)
C
      ELSE
C
CCC       WRITE(*,*) 'STMOVE: Resetting stagnation point'
C
C----- set new BL position -> panel position  pointers
       CALL IBLPAN(xfd)
C
C----- set new inviscid BL edge velocity UINV from QINV
       CALL UICALC(xfd)
C
C----- recalculate new arc length array
       CALL XICALC(xfd)
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS(xfd)
C
       IF(xfd%IST.GT.ISTOLD) THEN
C------ increase in number of points on top side (IS=1)
        IDIF = xfd%IST-ISTOLD
C
        xfd%ITRAN(1) = xfd%ITRAN(1) + IDIF
        xfd%ITRAN(2) = xfd%ITRAN(2) - IDIF
C
C------ move top side BL variables downstream
        DO 110 IBL=xfd%NBL(1), IDIF+2, -1
          xfd%CTAU(IBL,1) = xfd%CTAU(IBL-IDIF,1)
          xfd%THET(IBL,1) = xfd%THET(IBL-IDIF,1)
          xfd%DSTR(IBL,1) = xfd%DSTR(IBL-IDIF,1)
          xfd%UEDG(IBL,1) = xfd%UEDG(IBL-IDIF,1)
  110   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = xfd%UEDG(IDIF+2,1)/xfd%XSSI(IDIF+2,1)
        DO 115 IBL=IDIF+1, 2, -1
          xfd%CTAU(IBL,1) = xfd%CTAU(IDIF+2,1)
          xfd%THET(IBL,1) = xfd%THET(IDIF+2,1)
          xfd%DSTR(IBL,1) = xfd%DSTR(IDIF+2,1)
          xfd%UEDG(IBL,1) = DUDX * xfd%XSSI(IBL,1)
  115   CONTINUE
C
C------ move bottom side BL variables upstream
        DO 120 IBL=2, xfd%NBL(2)
          xfd%CTAU(IBL,2) = xfd%CTAU(IBL+IDIF,2)
          xfd%THET(IBL,2) = xfd%THET(IBL+IDIF,2)
          xfd%DSTR(IBL,2) = xfd%DSTR(IBL+IDIF,2)
          xfd%UEDG(IBL,2) = xfd%UEDG(IBL+IDIF,2)
  120   CONTINUE            
C
       ELSE
C------ increase in number of points on bottom side (IS=2)
        IDIF = ISTOLD-xfd%IST
C
        xfd%ITRAN(1) = xfd%ITRAN(1) - IDIF
        xfd%ITRAN(2) = xfd%ITRAN(2) + IDIF
C
C------ move bottom side BL variables downstream
        DO 210 IBL=xfd%NBL(2), IDIF+2, -1
          xfd%CTAU(IBL,2) = xfd%CTAU(IBL-IDIF,2)
          xfd%THET(IBL,2) = xfd%THET(IBL-IDIF,2)
          xfd%DSTR(IBL,2) = xfd%DSTR(IBL-IDIF,2)
          xfd%UEDG(IBL,2) = xfd%UEDG(IBL-IDIF,2)
  210   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = xfd%UEDG(IDIF+2,2)/xfd%XSSI(IDIF+2,2)


c        write(*,*) 'idif Ue xi dudx', 
c     &    idif, UEDG(idif+2,2), xssi(idif+2,2), dudx

        DO 215 IBL=IDIF+1, 2, -1
          xfd%CTAU(IBL,2) = xfd%CTAU(IDIF+2,2)
          xfd%THET(IBL,2) = xfd%THET(IDIF+2,2)
          xfd%DSTR(IBL,2) = xfd%DSTR(IDIF+2,2)
          xfd%UEDG(IBL,2) = DUDX * xfd%XSSI(IBL,2)
  215   CONTINUE

c        write(*,*) 'Uenew xinew', idif+1, uedg(idif+1,2), xssi(idif+1,2)

C
C------ move top side BL variables upstream
        DO 220 IBL=2, xfd%NBL(1)
          xfd%CTAU(IBL,1) = xfd%CTAU(IBL+IDIF,1)
          xfd%THET(IBL,1) = xfd%THET(IBL+IDIF,1)
          xfd%DSTR(IBL,1) = xfd%DSTR(IBL+IDIF,1)
          xfd%UEDG(IBL,1) = xfd%UEDG(IBL+IDIF,1)
  220   CONTINUE            
       ENDIF
C
C----- tweak Ue so it's not zero, in case stag. point is right on node
       UEPS = 1.0E-7
       DO IS = 1, 2
         DO IBL = 2, xfd%NBL(IS)
           I = xfd%IPAN(IBL,IS)
           IF(xfd%UEDG(IBL,IS).LE.UEPS) THEN
            xfd%UEDG(IBL,IS) = UEPS
            xfd%QVIS(I) = xfd%VTI(IBL,IS)*UEPS
            xfd%GAM(I)  = xfd%VTI(IBL,IS)*UEPS
           ENDIF
         ENDDO
       ENDDO
C
      ENDIF
C
C---- set new mass array since Ue has been tweaked
      DO 50 IS=1, 2
        DO 510 IBL=2, xfd%NBL(IS)
          xfd%MASS(IBL,IS) = xfd%DSTR(IBL,IS)*xfd%UEDG(IBL,IS)
  510   CONTINUE
   50 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Converges viscous operating point
C
C===================================================================70
      SUBROUTINE VISCAL(xfd,bld,xbd,NITER1)

      use xfoil_data_mod
      type(xfoil_data_type), intent(inout) :: xfd
      type(blpar_data_type), intent(inout) :: bld
      type(xbl_data_type), intent(inout) :: xbd
C
C---- convergence tolerance
      DATA EPS1 / 1.0E-4 /
C
      NITER = NITER1

C     DP mod: variable to notify of infinite loop condition (and halt)
      xfd%XFOIL_FAIL = .FALSE.
C
C---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.xfd%LWAKE) THEN
       CALL XYWAKE(xfd)
      ENDIF
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC(xfd)
C
C---- set velocities on airfoil and wake for initial alpha
      CALL QISET(xfd)
C
      IF(.NOT.xfd%LIPAN) THEN
C
       IF(xfd%LBLINI) CALL GAMQV(xfd)
C
C----- locate stagnation point arc length position and panel index
       CALL STFIND(xfd)
C
C----- set  BL position -> panel position  pointers
       CALL IBLPAN(xfd)
C
C----- calculate surface arc length array for current stagnation point location
       CALL XICALC(xfd)
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS(xfd)
C
      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC(xfd)
C
      IF(.NOT.xfd%LBLINI) THEN
C
C----- set initial Ue from inviscid Ue
       DO IBL=1, xfd%NBL(1)
         xfd%UEDG(IBL,1) = xfd%UINV(IBL,1)
       ENDDO
C
       DO IBL=1, xfd%NBL(2)
         xfd%UEDG(IBL,2) = xfd%UINV(IBL,2)
       ENDDO
C
      ENDIF
C
      IF(xfd%LVCONV) THEN
C----- set correct CL if converged point exists
       CALL QVFUE(xfd)
       IF(xfd%LVISC) THEN
        CALL CPCALC(xfd%N+xfd%NW,xfd%QVIS,xfd%QINF,xfd%MINF,xfd%CPV
     &  ,xfd%SILENT_MODE)
        CALL CPCALC(xfd%N+xfd%NW,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
       ELSE
        CALL CPCALC(xfd%N,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
       ENDIF
       CALL GAMQV(xfd)
       CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA,xfd%MINF
     &  ,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &             xfd%CL,xfd%CM,xfd%CDP, xfd%CL_ALF,xfd%CL_MSQ)
       CALL CDCALC(xfd)
      ENDIF
C
C---- set up source influence matrix if it doesn't exist
      IF(.NOT.xfd%LWDIJ .OR. .NOT.xfd%LADIJ) CALL QDCALC(xfd)
C
C---- Newton iteration for entire BL solution
C     DP mod: set default NITER to 10
C      IF(NITER.EQ.0) CALL ASKI('Enter number of iterations^',NITER)
      IF(NITER.EQ.0) NITER = 10
C     DP mod: added SILENT_MODE option
      IF (.NOT. xfd%SILENT_MODE) THEN
        WRITE(*,*)
        WRITE(*,*) 'Solving BL system ...'
      END IF
      DO 1000 ITER=1, NITER
C
C------ fill Newton system for BL variables
        CALL SETBL(xfd,bld,xbd)
C       DP mod: check for infinite loop condition
        IF (xfd%XFOIL_FAIL) THEN
          xfd%CL = -0.1
          xfd%CD = 1000.0
          xfd%CM = -10.0
          xfd%RMSBL = 1000.0
          RETURN
        ENDIF
C
C------ solve Newton system with custom solver
        CALL BLSOLV(xfd)
C
C------ update BL variables
        CALL UPDATE(xfd)
C
        IF(xfd%LALFA) THEN
C------- set new freestream Mach, Re from new CL
         CALL MRCL(xfd,xfd%CL,xfd%MINF_CL,REINF_CL)
         CALL COMSET(xfd)
        ELSE
C------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET(xfd)
         CALL UICALC(xfd)
        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE(xfd)
C
C------ set GAM distribution from QVIS
        CALL GAMQV(xfd)
C
C------ relocate stagnation point
        CALL STMOVE(xfd)
C
C------ set updated CL,CD
        CALL CLCALC(xfd%N,xfd%X,xfd%Y,xfd%GAM,xfd%GAM_A,xfd%ALFA
     &  ,xfd%MINF,xfd%QINF, xfd%XCMREF,xfd%YCMREF,
     &              xfd%CL,xfd%CM,xfd%CDP,xfd%CL_ALF,xfd%CL_MSQ)
        CALL CDCALC(xfd)
C
C------ display changes and test for convergence
C       DP mod: added SILENT_MODE options
        IF(.NOT. xfd%SILENT_MODE .AND. xfd%RLX.LT.1.0) 
     &   WRITE(*,2000) ITER, xfd%RMSBL, xfd%RMXBL, xfd%VMXBL,xfd%IMXBL
     &  ,xfd%ISMXBL,xfd%RLX
        IF(.NOT. xfd%SILENT_MODE .AND. xfd%RLX.EQ.1.0) 
     &   WRITE(*,2010) ITER, xfd%RMSBL, xfd%RMXBL, xfd%VMXBL,xfd%IMXBL
     &  ,xfd%ISMXBL
         CDPDIF = xfd%CD - xfd%CDF
         IF(.NOT. xfd%SILENT_MODE)
     &     WRITE(*,2020) xfd%ALFA/xfd%DTOR, xfd%CL, xfd%CM, xfd%CD,
     &   xfd%CDF, CDPDIF
C
        IF(xfd%RMSBL .LT. EPS1) THEN
         xfd%LVCONV = .TRUE.
         xfd%AVISC = xfd%ALFA
         xfd%MVISC = xfd%MINF
         GO TO 90
        ENDIF
C
 1000 CONTINUE
C     DP mod: added SILENT_MODE option
      IF(.NOT. xfd%SILENT_MODE) WRITE(*,*) 'VISCAL:  Convergence failed'
     &  
C
   90 CONTINUE
      CALL CPCALC(xfd%N+xfd%NW,xfd%QINV,xfd%QINF,xfd%MINF,xfd%CPI
     &  ,xfd%SILENT_MODE)
      CALL CPCALC(xfd%N+xfd%NW,xfd%QVIS,xfd%QINF,xfd%MINF,xfd%CPV
     &  ,xfd%SILENT_MODE)
C      IF(LFLAP) CALL MHINGE
      RETURN
C....................................................................
 2000   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,
     &     '   RLX:',F6.3)
 2010   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /
     &     1X,3X,'  Cm =', F8.4, '     CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL
