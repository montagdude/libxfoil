/* This file is part of libxfoil.
 
   libxfoil is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   libxfoil is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with libxfoil.  If not, see <http://www.gnu.org/licenses/>.
 
   Copyright (C) 2018 Daniel Prosser

   See xfoil_interface.f90 for descriptions of inputs and outputs. */

#pragma once

/* Xfoil dimensioning constants */

#define IQX 360
#define ISX 2
#define IBX 1440
#define IWX 47
#define IZX 407
#define IVX 277
#define NCOM 73

#include <stdbool.h>

/* Xfoil settings structs */

typedef struct
{
  double ncrit;
  double xtript, xtripb;
  bool viscous_mode;
  bool silent_mode;
  int maxit;
  double vaccel;
} xfoil_options_type;

typedef struct
{
  int npan;
  double cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2;
} xfoil_geom_options_type;

/* Data structs. Users only need to maintain an instance of the xfoil_data_group
   type, though the subtypes also need to be defined here. */

typedef struct
{
  double PI, HOPI, QOPI, DTOR;
  bool SILENT_MODE, VISCOUS_MODE;
  int MAXIT;
  double GAM[IQX], GAMU[IQX*2], QINVU[IZX*2], GAM_A[IQX];
  double QINV[IZX], QINV_A[IZX];
  bool LGAMU, LQAIJ, SHARP, LVISC, LWAKE, LVCONV, LWDIJ, LIPAN;
  bool LBLINI, LADIJ, LALFA;
  int RETYP, MATYP, ITMAX;
  double AIJ[IQX*IQX], BIJ[IQX*IZX], DIJ[IZX*IZX], CIJ[IWX*IQX];
  double DZDG[IQX], DZDN[IQX], DQDG[IQX], DZDM[IZX], DQDM[IZX];
  double X[IZX], Y[IZX], NX[IZX], NY[IZX], S[IZX], APANEL[IZX];
  double SIG[IZX], XP[IZX], YP[IZX];
  int N, NB, NPAN, NW, IST, NSYS;
  double PSIO, QINF, ALFA, Z_QINF, Z_ALFA, Z_QDOF0, Z_QDOF1;
  double Z_QDOF2, Z_QDOF3, ANTE, ASTE, DSTE, ADEG, AMAX;
  double QF0[IQX], QF1[IQX], QF2[IQX], QF3[IQX];
  int AIJPIV[IQX], IBLTE[ISX], NBL[ISX];
  int IPAN[IVX*ISX], ISYS[IVX*ISX];
  double SIGTE, GAMTE, SIGTE_A, GAMTE_A, MINF, MINF1, REINF, REINF1;
  double TKLAM, TKL_MSQ, CPSTAR, QSTAR, GAMMA, GAMM1;
  double XCMREF, YCMREF, CL, CM, CD, CDP, CDF, CL_ALF, CL_MSQ, SBLE;
  double XB[IBX], YB[IBX], SB[IBX], XBP[IBX], YBP[IBX], SNEW[5*IBX];
  double W1[6*IQX], W2[6*IQX], W3[6*IQX];
  double W4[6*IQX], W5[6*IQX], W6[6*IQX];
  double XLE, YLE, XTE, YTE, CHORD, SLE;
  double CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1, XPREF2;
  double MINF_CL, COSA, SINA, ACRIT, RLX, VACCEL;
  double CPI[IZX], CPV[IZX], QVIS[IZX];
  double VTI[IVX*ISX], XSSI[IVX*ISX];
  double AWAKE, AVISC, MVISC, CLSPEC, QTAN1, QTAN2, SST, SST_GO;
  double SST_GP;
  double WGAP[IWX], XSTRIP[ISX], XSSITR[ISX];
  double UINV[IVX*ISX], UINV_A[IVX*ISX], UEDG[IVX*ISX];
  double THET[IVX*ISX], DSTR[IVX*ISX], CTAU[IVX*ISX];
  double MASS[IVX*ISX], TAU[IVX*ISX], DIS[IVX*ISX], CTQ[IVX*ISX];
  double DELT[IVX*ISX], TSTR[IVX*ISX], USLP[IVX*ISX];
  int IDAMP;
  int ITRAN[ISX], IMXBL, ISMXBL;
  bool TFORCE[ISX];
  double VM[3*IZX*IZX], VA[3*2*IZX], VB[3*2*IZX], VDEL[3*2*IZX];
  double VZ[3*2], XOCTR[ISX], YOCTR[ISX];
  double RMSBL, RMXBL, WAKLEN;
  double UNEW[IVX*2], U_AC[IVX*2];
  double QNEW[IQX], Q_AC[IQX];
  char VMXBL;
  double THICKB, XTHICKB, THICKM, XTHICKM, CAMBR, XCAMBR;
  bool XFOIL_FAIL;
} xfoil_data_type;

typedef struct
{
  double SCCON, GACON, GBCON, GCCON, DLCON, CTRCON, CTRCEX, DUXCON;
  double CTCON, CFFAC;
} blpar_data_type;

typedef struct
{
  int IDAMPV;
  bool SIMI, TRAN, TURB, WAKE;
  bool TRFORC, TRFREE;
  double X1, U1, T1, D1, S1, AMPL1, U1_UEI, U1_MS, DW1;
  double H1, H1_T1, H1_D1;
  double M1, M1_U1, M1_MS;
  double R1, R1_U1, R1_MS;
  double V1, V1_U1, V1_MS, V1_RE;
  double HK1, HK1_U1, HK1_T1, HK1_D1, HK1_MS;
  double HS1, HS1_U1, HS1_T1, HS1_D1, HS1_MS, HS1_RE;
  double HC1, HC1_U1, HC1_T1, HC1_D1, HC1_MS;
  double RT1, RT1_U1, RT1_T1, RT1_MS, RT1_RE;
  double CF1, CF1_U1, CF1_T1, CF1_D1, CF1_MS, CF1_RE;
  double DI1, DI1_U1, DI1_T1, DI1_D1, DI1_S1, DI1_MS, DI1_RE;
  double US1, US1_U1, US1_T1, US1_D1, US1_MS, US1_RE;
  double CQ1, CQ1_U1, CQ1_T1, CQ1_D1, CQ1_MS, CQ1_RE;
  double DE1, DE1_U1, DE1_T1, DE1_D1, DE1_MS;
  double X2, U2, T2, D2, S2, AMPL2, U2_UEI, U2_MS, DW2;
  double H2, H2_T2, H2_D2;
  double M2, M2_U2, M2_MS;
  double R2, R2_U2, R2_MS;
  double V2, V2_U2, V2_MS, V2_RE;
  double HK2, HK2_U2, HK2_T2, HK2_D2, HK2_MS;
  double HS2, HS2_U2, HS2_T2, HS2_D2, HS2_MS, HS2_RE;
  double HC2, HC2_U2, HC2_T2, HC2_D2, HC2_MS;
  double RT2, RT2_U2, RT2_T2, RT2_MS, RT2_RE;
  double CF2, CF2_U2, CF2_T2, CF2_D2, CF2_MS, CF2_RE;
  double DI2, DI2_U2, DI2_T2, DI2_D2, DI2_S2, DI2_MS, DI2_RE;
  double US2, US2_U2, US2_T2, US2_D2, US2_MS, US2_RE;
  double CQ2, CQ2_U2, CQ2_T2, CQ2_D2, CQ2_MS, CQ2_RE;
  double DE2, DE2_U2, DE2_T2, DE2_D2, DE2_MS;
  double CFM, CFM_MS, CFM_RE;
  double CFM_U1, CFM_T1, CFM_D1, CFM_U2, CFM_T2, CFM_D2;
  double XT, XT_A1, XT_MS, XT_RE, XT_XF;
  double XT_X1, XT_T1, XT_D1, XT_U1, XT_X2, XT_T2, XT_D2, XT_U2;
  double C1SAV[NCOM], C2SAV[NCOM];
  double DWTE, QINFBL, TKBL, TKBL_MS, RSTBL, RSTBL_MS, HSTINV;
  double HSTINV_MS, REYBL, REYBL_MS, REYBL_RE, GAMBL, GM1BL;
  double HVRA, BULE, XIFORC, AMCRIT;
  double VS1[4*5], VS2[4*5];
  double VSREZ[4], VSR[4], VSM[4], VSX[4];
} xbl_data_type;

typedef struct
{
  xfoil_data_type xfd;
  blpar_data_type bld;
  xbl_data_type xbd;
} xfoil_data_group;

extern void xfoil_defaults(xfoil_data_group *xdg,
                           const xfoil_options_type *xfoil_options);
extern void xfoil_set_paneling(xfoil_data_group *xdg,
                               const xfoil_geom_options_type *geom_opts);
extern void xfoil_set_airfoil(xfoil_data_group *xdg, const double xin[],
                              const double zin[], const int *npointin);
extern void xfoil_smooth_paneling(xfoil_data_group *xdg, int *stat);
extern void xfoil_apply_flap_deflection(xfoil_data_group *xdg,
                                        const double *xflap,
                                        const double *yflap,
                                        const int *y_flap_spec,
                                        const double *degrees, int *npointout,
                                        int *stat);
extern void xfoil_modify_tegap(xfoil_data_group *xdg, const double *gap,
                               const double *blendloc, int *npointout,
                               int *stat);
extern void xfoil_get_airfoil(const xfoil_data_group *xdg, double xout[],
                              double zout[], const int *npoint, int *stat);
extern void xfoil_geometry_info(const xfoil_data_group *xdg, double *maxt,
                                double *xmaxt, double *maxc, double *xmaxc,
                                int *stat);
extern void xfoil_set_reynolds_number(xfoil_data_group *xdg, const double *re);
extern void xfoil_set_mach_number(xfoil_data_group *xdg, const double *mach);
extern void xfoil_reinitialize_bl(xfoil_data_group *xdg);
extern void xfoil_specal(xfoil_data_group *xdg, const double *alpha_spec,
                         double *alpha, double *lift, double *drag,
                         double *moment, bool *converged, int *stat);
extern void xfoil_speccl(xfoil_data_group *xdg, const double *cl_spec,
                         double *alpha, double *lift, double *drag,
                         double *moment, bool *converged, int *stat);
extern void xfoil_get_transloc(const xfoil_data_group *xdg, double *xtranst,
                               double *ztranst, double *xtransb,
                               double *ztransb);
extern void xfoil_get_cp(const xfoil_data_group *xdg, const int *npoint,
                         double cp[]);
extern void xfoil_get_cf(const xfoil_data_group *xdg, const int *npoint,
                         double cf[]);
extern void xfoil_get_uedge(const xfoil_data_group *xdg, const int *npoint,
                            double uedge[]);
extern void xfoil_get_deltastar(const xfoil_data_group *xdg, const int *npoint,
                                double deltastar[]);
extern void xfoil_get_diss(const xfoil_data_group *xdg, const int *npoint,
                           double diss[]);
extern void xfoil_get_hk(const xfoil_data_group *xdg, const int *npoint,
                         double hk[]);
extern void xfoil_get_retheta(const xfoil_data_group *xdg, const int *npoint,
                              double retheta[]);

/* The following methods utilize xfoil functionality, performing some
   calculations and returning a result, without needing to maintain an
   xfoil_data_group instance. */
extern void run_xfoil(const int *npointin, const double xin[],
                      const double zin[],
                      const xfoil_geom_options_type *geom_opts,
                      const int *noppoint, const double operating_points[],
                      const int op_modes[], const double reynolds_numbers[],
                      const double mach_numbers[], const bool *use_flap,
                      const double *x_flap, const double *y_flap,
                      const int *y_flap_spec, const double flap_degrees[],
                      const xfoil_options_type *xfoil_opts, 
                      const bool *reinitialize, const bool *fix_unconverged,
                      double lift[], double drag[], double moment[],
                      double viscrms[], double alpha[], double xtrt[],
                      double xtrb[], int *stat, const double ncrit_per_point[]);
extern void naca_4_digit(const char des[4], const int *npointside,
                         double xout[], double zout[], int *nout);
extern void naca_5_digit(const char des[5], const int *npointside,
                         double xout[], double zout[], int *nout, int *stat);
extern void xfoil_spline_coordinates(const double x[], const double z[],
                                     const int *npt, double s[], double xs[],
                                     double zs[]);
extern void xfoil_eval_spline(const double x[], const double z[],
                              const double s[], const double xs[],
                              const double zs[], const int *npt,
                              const double *sc, double *xc, double *zc);
extern void xfoil_lefind(const double x[], const double z[], const double s[],
                         const double xs[], const double zs[], const int *npt,
                         double *sle, double *xle, double *zle);
