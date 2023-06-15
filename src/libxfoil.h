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

   Copyright (C) 2019 Daniel Prosser

   See libxfoil.f90 for descriptions of inputs and outputs. */

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
  double *GAM, *GAMU, *QINVU, *GAM_A;
  double *QINV, *QINV_A;
  bool LGAMU, LQAIJ, SHARP, LVISC, LWAKE, LVCONV, LWDIJ, LIPAN;
  bool LBLINI, LADIJ, LALFA;
  int RETYP, MATYP, ITMAX;
  double *AIJ, *BIJ, *DIJ, *CIJ;
  double *DZDG, *DZDN, *DQDG, *DZDM, *DQDM;
  double *X, *Y, *NX, *NY, *S, *APANEL;
  double *SIG, *XP, *YP;
  int N, NB, NPAN, NW, IST, NSYS;
  double PSIO, QINF, ALFA, Z_QINF, Z_ALFA, Z_QDOF0, Z_QDOF1;
  double Z_QDOF2, Z_QDOF3, ANTE, ASTE, DSTE, ADEG, AMAX;
  double *QF0, *QF1, *QF2, *QF3;
  int *AIJPIV, *IBLTE, *NBL;
  int *IPAN, *ISYS;
  double SIGTE, GAMTE, SIGTE_A, GAMTE_A, MINF, MINF1, REINF, REINF1;
  double TKLAM, TKL_MSQ, CPSTAR, QSTAR, GAMMA, GAMM1;
  double XCMREF, YCMREF, CL, CM, CD, CDP, CDF, CL_ALF, CL_MSQ, SBLE;
  double *XB, *YB, *SB, *XBP, *YBP, *SNEW;
  double *W1, *W2, *W3;
  double *W4, *W5, *W6;
  double XLE, YLE, XTE, YTE, CHORD, SLE;
  double CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1, XPREF2;
  double MINF_CL, COSA, SINA, ACRIT, RLX, VACCEL;
  double *CPI, *CPV, *QVIS;
  double *VTI, *XSSI;
  double AWAKE, AVISC, MVISC, CLSPEC, QTAN1, QTAN2, SST, SST_GO;
  double SST_GP;
  double *WGAP, *XSTRIP, *XSSITR;
  double *UINV, *UINV_A, *UEDG;
  double *THET, *DSTR, *CTAU;
  double *MASS, *TAU, *DIS, *CTQ;
  double *DELT, *TSTR, *USLP;
  int IDAMP;
  int *ITRAN;
  int IMXBL, ISMXBL;
  bool *TFORCE;
  double *VM, *VA, *VB, *VDEL;
  double *VZ, *XOCTR, *YOCTR;
  double RMSBL, RMXBL, WAKLEN;
  double *UNEW, *U_AC;
  double *QNEW, *Q_AC;
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
  double *C1SAV, *C2SAV;
  double DWTE, QINFBL, TKBL, TKBL_MS, RSTBL, RSTBL_MS, HSTINV;
  double HSTINV_MS, REYBL, REYBL_MS, REYBL_RE, GAMBL, GM1BL;
  double HVRA, BULE, XIFORC, AMCRIT;
  double *VS1, *VS2;
  double *VSREZ, *VSR, *VSM, *VSX;
} xbl_data_type;

typedef struct
{
  xfoil_data_type xfd;
  blpar_data_type bld;
  xbl_data_type xbd;
} xfoil_data_group;

/* Memory management */

extern void xfoil_init(xfoil_data_group *xdg);
extern void xfoil_cleanup(xfoil_data_group *xdg);
extern void xfoil_copy(const xfoil_data_group *xdg_from,
                       xfoil_data_group *xdg_to);

/* Xfoil routines modifying xfoil_data_group */

extern void xfoil_defaults(xfoil_data_group *xdg,
                           const xfoil_options_type *xfoil_options);
extern void xfoil_set_paneling(xfoil_data_group *xdg,
                               const xfoil_geom_options_type *geom_opts);
extern void xfoil_set_buffer_airfoil(xfoil_data_group *xdg, const double xin[],
                                     const double zin[], const int *npointin);
extern void xfoil_get_buffer_airfoil(const xfoil_data_group *xdg, double xout[],
                                     double zout[], const int *npoint,
                                     int *stat);
extern void xfoil_get_current_airfoil(const xfoil_data_group *xdg,
                                      double xout[], double zout[],
                                      const int *npoint, int *stat);
extern void xfoil_smooth_paneling(xfoil_data_group *xdg, int *stat);
extern void xfoil_update_current_airfoil(const xfoil_data_group *xdg, int *stat);
extern void xfoil_apply_flap_deflection(xfoil_data_group *xdg,
                                        const double *xflap,
                                        const double *zflap,
                                        const int *z_flap_spec,
                                        const double *degrees, int *npointout,
                                        int *stat);
extern void xfoil_modify_tegap(xfoil_data_group *xdg, const double *gap,
                               const double *blendloc, int *npointout,
                               int *stat);
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
extern void xfoil_get_ampl(const xfoil_data_group *xdg, const int *npoint,
                           double ampl[]);

/* Wake data (also uses xfoil_data_group) */
extern void xfoil_get_wakepoints(const xfoil_data_group *xdg, int *nwake);
extern void xfoil_get_wake_geometry(const xfoil_data_group *xdg,
                                    const int *nwake, double xw[],
                                    double zw[]);
extern void xfoil_get_wake_cp(const xfoil_data_group *xdg, const int *nwake,
                              double cp[]);
extern void xfoil_get_wake_uedge(const xfoil_data_group *xdg, const int *nwake,
                                 double uedge[]);
extern void xfoil_get_wake_deltastar(const xfoil_data_group *xdg,
                                     const int *nwake, double deltastar[]);

/* The following methods utilize xfoil functionality, performing some
   calculations and returning a result, without needing to maintain an
   xfoil_data_group instance. */
extern void run_xfoil(const int *npointin, const double xin[],
                      const double zin[],
                      const xfoil_geom_options_type *geom_opts,
                      const int *noppoint, const double operating_points[],
                      const int op_modes[], const double reynolds_numbers[],
                      const double mach_numbers[], const bool *use_flap,
                      const double *x_flap, const double *z_flap,
                      const int *z_flap_spec, const double flap_degrees[],
                      const xfoil_options_type *xfoil_opts,
                      const bool *reinitialize, const bool *fix_unconverged,
                      double lift[], double drag[], double moment[],
                      double viscrms[], double alpha[], double xtrt[],
                      double xtrb[], int *stat);
extern void naca_4_digit(const double *camber, const double *xcamber,
                         const double *thick, const int *npointside,
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
