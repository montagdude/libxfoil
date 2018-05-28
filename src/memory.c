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

#include <stdbool.h>
#include <stdlib.h>
#include "xfoil_interface.h"

/******************************************************************************/
//
// Allocates memory for structs in xfoil_data_group. This is the C backend
// for the xfoil_init library method.
//
/******************************************************************************/
void allocate_xdg(xfoil_data_group *xdg)
{
  xdg->xfd.GAM = malloc(IQX*sizeof(double));
  xdg->xfd.GAMU = malloc(IQX*2*sizeof(double));
  xdg->xfd.QINVU = malloc(IZX*2*sizeof(double));
  xdg->xfd.GAM_A = malloc(IQX*sizeof(double));
  xdg->xfd.QINV = malloc(IZX*sizeof(double));
  xdg->xfd.QINV_A = malloc(IZX*sizeof(double));
  xdg->xfd.AIJ = malloc(IQX*IQX*sizeof(double));
  xdg->xfd.BIJ = malloc(IQX*IZX*sizeof(double));
  xdg->xfd.DIJ = malloc(IZX*IZX*sizeof(double));
  xdg->xfd.CIJ = malloc(IWX*IQX*sizeof(double));
  xdg->xfd.DZDG = malloc(IQX*sizeof(double));
  xdg->xfd.DZDN = malloc(IQX*sizeof(double));
  xdg->xfd.DQDG = malloc(IQX*sizeof(double));
  xdg->xfd.DZDM = malloc(IZX*sizeof(double));
  xdg->xfd.DQDM = malloc(IZX*sizeof(double));
  xdg->xfd.X = malloc(IZX*sizeof(double));
  xdg->xfd.Y = malloc(IZX*sizeof(double));
  xdg->xfd.NX = malloc(IZX*sizeof(double));
  xdg->xfd.NY = malloc(IZX*sizeof(double));
  xdg->xfd.S = malloc(IZX*sizeof(double));
  xdg->xfd.APANEL = malloc(IZX*sizeof(double));
  xdg->xfd.SIG = malloc(IZX*sizeof(double));
  xdg->xfd.XP = malloc(IZX*sizeof(double));
  xdg->xfd.YP = malloc(IZX*sizeof(double));
  xdg->xfd.QF0 = malloc(IQX*sizeof(double));
  xdg->xfd.QF1 = malloc(IQX*sizeof(double));
  xdg->xfd.QF2 = malloc(IQX*sizeof(double));
  xdg->xfd.QF3 = malloc(IQX*sizeof(double));
  xdg->xfd.AIJPIV = malloc(IQX*sizeof(int));
  xdg->xfd.IBLTE = malloc(ISX*sizeof(int));
  xdg->xfd.NBL = malloc(ISX*sizeof(int));
  xdg->xfd.IPAN = malloc(IVX*ISX*sizeof(int));
  xdg->xfd.ISYS = malloc(IVX*ISX*sizeof(int));
  xdg->xfd.XB = malloc(IBX*sizeof(double));
  xdg->xfd.YB = malloc(IBX*sizeof(double));
  xdg->xfd.SB = malloc(IBX*sizeof(double));
  xdg->xfd.XBP = malloc(IBX*sizeof(double));
  xdg->xfd.YBP = malloc(IBX*sizeof(double));
  xdg->xfd.SNEW = malloc(5*IBX*sizeof(double));
  xdg->xfd.W1 = malloc(6*IQX*sizeof(double));
  xdg->xfd.W2 = malloc(6*IQX*sizeof(double));
  xdg->xfd.W3 = malloc(6*IQX*sizeof(double));
  xdg->xfd.W4 = malloc(6*IQX*sizeof(double));
  xdg->xfd.W5 = malloc(6*IQX*sizeof(double));
  xdg->xfd.W6 = malloc(6*IQX*sizeof(double));
  xdg->xfd.CPI = malloc(IZX*sizeof(double));
  xdg->xfd.CPV = malloc(IZX*sizeof(double));
  xdg->xfd.QVIS = malloc(IZX*sizeof(double));
  xdg->xfd.VTI = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.XSSI = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.WGAP = malloc(IWX*sizeof(double));
  xdg->xfd.XSTRIP = malloc(ISX*sizeof(double));
  xdg->xfd.XSSITR = malloc(ISX*sizeof(double));
  xdg->xfd.UINV = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.UINV_A = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.UEDG = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.THET = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.DSTR = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.CTAU = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.MASS = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.TAU = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.DIS = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.CTQ = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.DELT = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.TSTR = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.USLP = malloc(IVX*ISX*sizeof(double));
  xdg->xfd.ITRAN = malloc(ISX*sizeof(int));
  xdg->xfd.TFORCE = malloc(ISX*sizeof(bool));
  xdg->xfd.VM = malloc(3*IZX*IZX*sizeof(double));
  xdg->xfd.VA = malloc(3*2*IZX*sizeof(double));
  xdg->xfd.VB = malloc(3*2*IZX*sizeof(double));
  xdg->xfd.VDEL = malloc(3*2*IZX*sizeof(double));
  xdg->xfd.VZ = malloc(3*2*sizeof(double));
  xdg->xfd.XOCTR = malloc(ISX*sizeof(double));
  xdg->xfd.YOCTR = malloc(ISX*sizeof(double));
  xdg->xfd.UNEW = malloc(IVX*2*sizeof(double));
  xdg->xfd.U_AC = malloc(IVX*2*sizeof(double));
  xdg->xfd.QNEW = malloc(IQX*sizeof(double));
  xdg->xfd.Q_AC = malloc(IQX*sizeof(double));
  xdg->xbd.C1SAV = malloc(NCOM*sizeof(double));
  xdg->xbd.C2SAV = malloc(NCOM*sizeof(double));
  xdg->xbd.VS1 = malloc(4*5*sizeof(double));
  xdg->xbd.VS2 = malloc(4*5*sizeof(double));
  xdg->xbd.VSREZ = malloc(4*sizeof(double));
  xdg->xbd.VSR = malloc(4*sizeof(double));
  xdg->xbd.VSM = malloc(4*sizeof(double));
  xdg->xbd.VSX = malloc(4*sizeof(double));
}

/******************************************************************************/
//
// Frees memory for structs in xfoil_data_group. This is the C backend for the
// for the xfoil_cleanup library method.
//
/******************************************************************************/
void free_xdg(xfoil_data_group *xdg)
{
  free(xdg->xfd.GAM);
  free(xdg->xfd.GAMU);
  free(xdg->xfd.QINVU);
  free(xdg->xfd.GAM_A);
  free(xdg->xfd.QINV);
  free(xdg->xfd.QINV_A);
  free(xdg->xfd.AIJ);
  free(xdg->xfd.BIJ);
  free(xdg->xfd.DIJ);
  free(xdg->xfd.CIJ);
  free(xdg->xfd.DZDG);
  free(xdg->xfd.DZDN);
  free(xdg->xfd.DQDG);
  free(xdg->xfd.DZDM);
  free(xdg->xfd.DQDM);
  free(xdg->xfd.X);
  free(xdg->xfd.Y);
  free(xdg->xfd.NX);
  free(xdg->xfd.NY);
  free(xdg->xfd.S);
  free(xdg->xfd.APANEL);
  free(xdg->xfd.SIG);
  free(xdg->xfd.XP);
  free(xdg->xfd.YP);
  free(xdg->xfd.QF0);
  free(xdg->xfd.QF1);
  free(xdg->xfd.QF2);
  free(xdg->xfd.QF3);
  free(xdg->xfd.AIJPIV);
  free(xdg->xfd.IBLTE);
  free(xdg->xfd.NBL);
  free(xdg->xfd.IPAN);
  free(xdg->xfd.ISYS);
  free(xdg->xfd.XB);
  free(xdg->xfd.YB);
  free(xdg->xfd.SB);
  free(xdg->xfd.XBP);
  free(xdg->xfd.YBP);
  free(xdg->xfd.SNEW);
  free(xdg->xfd.W1);
  free(xdg->xfd.W2);
  free(xdg->xfd.W3);
  free(xdg->xfd.W4);
  free(xdg->xfd.W5);
  free(xdg->xfd.W6);
  free(xdg->xfd.CPI);
  free(xdg->xfd.CPV);
  free(xdg->xfd.QVIS);
  free(xdg->xfd.VTI);
  free(xdg->xfd.XSSI);
  free(xdg->xfd.WGAP);
  free(xdg->xfd.XSTRIP);
  free(xdg->xfd.XSSITR);
  free(xdg->xfd.UINV);
  free(xdg->xfd.UINV_A);
  free(xdg->xfd.UEDG);
  free(xdg->xfd.THET);
  free(xdg->xfd.DSTR);
  free(xdg->xfd.CTAU);
  free(xdg->xfd.MASS);
  free(xdg->xfd.TAU);
  free(xdg->xfd.DIS);
  free(xdg->xfd.CTQ);
  free(xdg->xfd.DELT);
  free(xdg->xfd.TSTR);
  free(xdg->xfd.USLP);
  free(xdg->xfd.ITRAN);
  free(xdg->xfd.TFORCE);
  free(xdg->xfd.VM);
  free(xdg->xfd.VA);
  free(xdg->xfd.VB);
  free(xdg->xfd.VDEL);
  free(xdg->xfd.VZ);
  free(xdg->xfd.XOCTR);
  free(xdg->xfd.YOCTR);
  free(xdg->xfd.UNEW);
  free(xdg->xfd.U_AC);
  free(xdg->xfd.QNEW);
  free(xdg->xfd.Q_AC);
  free(xdg->xbd.C1SAV);
  free(xdg->xbd.C2SAV);
  free(xdg->xbd.VS1);
  free(xdg->xbd.VS2);
  free(xdg->xbd.VSREZ);
  free(xdg->xbd.VSR);
  free(xdg->xbd.VSM);
  free(xdg->xbd.VSX);
}

/******************************************************************************/
//
// Copies everything in xdg_from to xdg_to. xfoil_init must have been already
// called for both.
//
/******************************************************************************/
void copy_xdg(xfoil_data_group *xdg_from, xfoil_data_group *xdg_to)
{
  int i;

  xdg_to->xfd.PI = xdg_from->xfd.PI;
  xdg_to->xfd.HOPI = xdg_from->xfd.HOPI;
  xdg_to->xfd.QOPI = xdg_from->xfd.QOPI;
  xdg_to->xfd.DTOR = xdg_from->xfd.DTOR;
  xdg_to->xfd.SILENT_MODE = xdg_from->xfd.SILENT_MODE;
  xdg_to->xfd.VISCOUS_MODE = xdg_from->xfd.VISCOUS_MODE;
  xdg_to->xfd.MAXIT = xdg_from->xfd.MAXIT;
  xdg_to->xfd.LGAMU = xdg_from->xfd.LGAMU;
  xdg_to->xfd.LQAIJ = xdg_from->xfd.LQAIJ;
  xdg_to->xfd.SHARP = xdg_from->xfd.SHARP;
  xdg_to->xfd.LVISC = xdg_from->xfd.LVISC;
  xdg_to->xfd.LWAKE = xdg_from->xfd.LWAKE;
  xdg_to->xfd.LVCONV = xdg_from->xfd.LVCONV;
  xdg_to->xfd.LWDIJ = xdg_from->xfd.LWDIJ;
  xdg_to->xfd.LIPAN = xdg_from->xfd.LIPAN;
  xdg_to->xfd.LBLINI = xdg_from->xfd.LBLINI;
  xdg_to->xfd.LADIJ = xdg_from->xfd.LADIJ;
  xdg_to->xfd.LALFA = xdg_from->xfd.LALFA;
  xdg_to->xfd.RETYP = xdg_from->xfd.RETYP;
  xdg_to->xfd.MATYP = xdg_from->xfd.MATYP;
  xdg_to->xfd.ITMAX = xdg_from->xfd.ITMAX;
  xdg_to->xfd.N = xdg_from->xfd.N;
  xdg_to->xfd.NB = xdg_from->xfd.NB;
  xdg_to->xfd.NPAN = xdg_from->xfd.NPAN;
  xdg_to->xfd.NW = xdg_from->xfd.NW;
  xdg_to->xfd.IST = xdg_from->xfd.IST;
  xdg_to->xfd.NSYS = xdg_from->xfd.NSYS;
  xdg_to->xfd.PSIO = xdg_from->xfd.PSIO;
  xdg_to->xfd.QINF = xdg_from->xfd.QINF;
  xdg_to->xfd.ALFA = xdg_from->xfd.ALFA;
  xdg_to->xfd.Z_QINF = xdg_from->xfd.Z_QINF;
  xdg_to->xfd.Z_ALFA = xdg_from->xfd.Z_ALFA;
  xdg_to->xfd.Z_QDOF0 = xdg_from->xfd.Z_QDOF0;
  xdg_to->xfd.Z_QDOF1 = xdg_from->xfd.Z_QDOF1;
  xdg_to->xfd.Z_QDOF2 = xdg_from->xfd.Z_QDOF2;
  xdg_to->xfd.Z_QDOF3 = xdg_from->xfd.Z_QDOF3;
  xdg_to->xfd.ANTE = xdg_from->xfd.ANTE;
  xdg_to->xfd.ASTE = xdg_from->xfd.ASTE;
  xdg_to->xfd.DSTE = xdg_from->xfd.DSTE;
  xdg_to->xfd.ADEG = xdg_from->xfd.ADEG;
  xdg_to->xfd.AMAX = xdg_from->xfd.AMAX;
  xdg_to->xfd.SIGTE = xdg_from->xfd.SIGTE;
  xdg_to->xfd.GAMTE = xdg_from->xfd.GAMTE;
  xdg_to->xfd.SIGTE_A = xdg_from->xfd.SIGTE_A;
  xdg_to->xfd.GAMTE_A = xdg_from->xfd.GAMTE_A;
  xdg_to->xfd.MINF = xdg_from->xfd.MINF;
  xdg_to->xfd.MINF1 = xdg_from->xfd.MINF1;
  xdg_to->xfd.REINF = xdg_from->xfd.REINF;
  xdg_to->xfd.REINF1 = xdg_from->xfd.REINF1;
  xdg_to->xfd.TKLAM = xdg_from->xfd.TKLAM;
  xdg_to->xfd.TKL_MSQ = xdg_from->xfd.TKL_MSQ;
  xdg_to->xfd.CPSTAR = xdg_from->xfd.CPSTAR;
  xdg_to->xfd.QSTAR = xdg_from->xfd.QSTAR;
  xdg_to->xfd.GAMMA = xdg_from->xfd.GAMMA;
  xdg_to->xfd.GAMM1 = xdg_from->xfd.GAMM1;
  xdg_to->xfd.XCMREF = xdg_from->xfd.XCMREF;
  xdg_to->xfd.YCMREF = xdg_from->xfd.YCMREF;
  xdg_to->xfd.CL = xdg_from->xfd.CL;
  xdg_to->xfd.CM = xdg_from->xfd.CM;
  xdg_to->xfd.CD = xdg_from->xfd.CD;
  xdg_to->xfd.CDP = xdg_from->xfd.CDP;
  xdg_to->xfd.CDF = xdg_from->xfd.CDF;
  xdg_to->xfd.CL_ALF = xdg_from->xfd.CL_ALF;
  xdg_to->xfd.CL_MSQ = xdg_from->xfd.CL_MSQ;
  xdg_to->xfd.SBLE = xdg_from->xfd.SBLE;
  xdg_to->xfd.XLE = xdg_from->xfd.XLE;
  xdg_to->xfd.YLE = xdg_from->xfd.YLE;
  xdg_to->xfd.XTE = xdg_from->xfd.XTE;
  xdg_to->xfd.YTE = xdg_from->xfd.YTE;
  xdg_to->xfd.CHORD = xdg_from->xfd.CHORD;
  xdg_to->xfd.SLE = xdg_from->xfd.SLE;
  xdg_to->xfd.CVPAR = xdg_from->xfd.CVPAR;
  xdg_to->xfd.CTERAT = xdg_from->xfd.CTERAT;
  xdg_to->xfd.CTRRAT = xdg_from->xfd.CTRRAT;
  xdg_to->xfd.XSREF1 = xdg_from->xfd.XSREF1;
  xdg_to->xfd.XSREF2 = xdg_from->xfd.XSREF2;
  xdg_to->xfd.XPREF1 = xdg_from->xfd.XPREF1;
  xdg_to->xfd.XPREF2 = xdg_from->xfd.XPREF2;
  xdg_to->xfd.MINF_CL = xdg_from->xfd.MINF_CL;
  xdg_to->xfd.COSA = xdg_from->xfd.COSA;
  xdg_to->xfd.SINA = xdg_from->xfd.SINA;
  xdg_to->xfd.ACRIT = xdg_from->xfd.ACRIT;
  xdg_to->xfd.RLX = xdg_from->xfd.RLX;
  xdg_to->xfd.VACCEL = xdg_from->xfd.VACCEL;
  xdg_to->xfd.AWAKE = xdg_from->xfd.AWAKE;
  xdg_to->xfd.AVISC = xdg_from->xfd.AVISC;
  xdg_to->xfd.MVISC = xdg_from->xfd.MVISC;
  xdg_to->xfd.CLSPEC = xdg_from->xfd.CLSPEC;
  xdg_to->xfd.QTAN1 = xdg_from->xfd.QTAN1;
  xdg_to->xfd.QTAN2 = xdg_from->xfd.QTAN2;
  xdg_to->xfd.SST = xdg_from->xfd.SST;
  xdg_to->xfd.SST_GO = xdg_from->xfd.SST_GO;
  xdg_to->xfd.SST_GP = xdg_from->xfd.SST_GP;
  xdg_to->xfd.IDAMP = xdg_from->xfd.IDAMP;
  xdg_to->xfd.IMXBL = xdg_from->xfd.IMXBL;
  xdg_to->xfd.ISMXBL = xdg_from->xfd.ISMXBL;
  xdg_to->xfd.RMSBL = xdg_from->xfd.RMSBL;
  xdg_to->xfd.RMXBL = xdg_from->xfd.RMXBL;
  xdg_to->xfd.WAKLEN = xdg_from->xfd.WAKLEN;
  xdg_to->xfd.VMXBL = xdg_from->xfd.VMXBL;
  xdg_to->xfd.THICKB = xdg_from->xfd.THICKB;
  xdg_to->xfd.XTHICKB = xdg_from->xfd.XTHICKB;
  xdg_to->xfd.THICKM = xdg_from->xfd.THICKM;
  xdg_to->xfd.XTHICKM = xdg_from->xfd.XTHICKM;
  xdg_to->xfd.CAMBR = xdg_from->xfd.CAMBR;
  xdg_to->xfd.XCAMBR = xdg_from->xfd.XCAMBR;
  xdg_to->xfd.XFOIL_FAIL = xdg_from->xfd.XFOIL_FAIL;
  xdg_to->xbd.IDAMPV = xdg_from->xbd.IDAMPV;
  xdg_to->xbd.SIMI = xdg_from->xbd.SIMI;
  xdg_to->xbd.TRAN = xdg_from->xbd.TRAN;
  xdg_to->xbd.TURB = xdg_from->xbd.TURB;
  xdg_to->xbd.WAKE = xdg_from->xbd.WAKE;
  xdg_to->xbd.TRFORC = xdg_from->xbd.TRFORC;
  xdg_to->xbd.TRFREE = xdg_from->xbd.TRFREE;
  xdg_to->xbd.X1 = xdg_from->xbd.X1;
  xdg_to->xbd.U1 = xdg_from->xbd.U1;
  xdg_to->xbd.T1 = xdg_from->xbd.T1;
  xdg_to->xbd.D1 = xdg_from->xbd.D1;
  xdg_to->xbd.S1 = xdg_from->xbd.S1;
  xdg_to->xbd.AMPL1 = xdg_from->xbd.AMPL1;
  xdg_to->xbd.U1_UEI = xdg_from->xbd.U1_UEI;
  xdg_to->xbd.U1_MS = xdg_from->xbd.U1_MS;
  xdg_to->xbd.DW1 = xdg_from->xbd.DW1;
  xdg_to->xbd.H1 = xdg_from->xbd.H1;
  xdg_to->xbd.H1_T1 = xdg_from->xbd.H1_T1;
  xdg_to->xbd.H1_D1 = xdg_from->xbd.H1_D1;
  xdg_to->xbd.M1 = xdg_from->xbd.M1;
  xdg_to->xbd.M1_U1 = xdg_from->xbd.M1_U1;
  xdg_to->xbd.M1_MS = xdg_from->xbd.M1_MS;
  xdg_to->xbd.R1 = xdg_from->xbd.R1;
  xdg_to->xbd.R1_U1 = xdg_from->xbd.R1_U1;
  xdg_to->xbd.R1_MS = xdg_from->xbd.R1_MS;
  xdg_to->xbd.V1 = xdg_from->xbd.V1;
  xdg_to->xbd.V1_U1 = xdg_from->xbd.V1_U1;
  xdg_to->xbd.V1_MS = xdg_from->xbd.V1_MS;
  xdg_to->xbd.V1_RE = xdg_from->xbd.V1_RE;
  xdg_to->xbd.HK1 = xdg_from->xbd.HK1;
  xdg_to->xbd.HK1_U1 = xdg_from->xbd.HK1_U1;
  xdg_to->xbd.HK1_T1 = xdg_from->xbd.HK1_T1;
  xdg_to->xbd.HK1_D1 = xdg_from->xbd.HK1_D1;
  xdg_to->xbd.HK1_MS = xdg_from->xbd.HK1_MS;
  xdg_to->xbd.HS1 = xdg_from->xbd.HS1;
  xdg_to->xbd.HS1_U1 = xdg_from->xbd.HS1_U1;
  xdg_to->xbd.HS1_T1 = xdg_from->xbd.HS1_T1;
  xdg_to->xbd.HS1_D1 = xdg_from->xbd.HS1_D1;
  xdg_to->xbd.HS1_MS = xdg_from->xbd.HS1_MS;
  xdg_to->xbd.HS1_RE = xdg_from->xbd.HS1_RE;
  xdg_to->xbd.HC1 = xdg_from->xbd.HC1;
  xdg_to->xbd.HC1_U1 = xdg_from->xbd.HC1_U1;
  xdg_to->xbd.HC1_T1 = xdg_from->xbd.HC1_T1;
  xdg_to->xbd.HC1_D1 = xdg_from->xbd.HC1_D1;
  xdg_to->xbd.HC1_MS = xdg_from->xbd.HC1_MS;
  xdg_to->xbd.RT1 = xdg_from->xbd.RT1;
  xdg_to->xbd.RT1_U1 = xdg_from->xbd.RT1_U1;
  xdg_to->xbd.RT1_T1 = xdg_from->xbd.RT1_T1;
  xdg_to->xbd.RT1_MS = xdg_from->xbd.RT1_MS;
  xdg_to->xbd.RT1_RE = xdg_from->xbd.RT1_RE;
  xdg_to->xbd.CF1 = xdg_from->xbd.CF1;
  xdg_to->xbd.CF1_U1 = xdg_from->xbd.CF1_U1;
  xdg_to->xbd.CF1_T1 = xdg_from->xbd.CF1_T1;
  xdg_to->xbd.CF1_D1 = xdg_from->xbd.CF1_D1;
  xdg_to->xbd.CF1_MS = xdg_from->xbd.CF1_MS;
  xdg_to->xbd.CF1_RE = xdg_from->xbd.CF1_RE;
  xdg_to->xbd.DI1 = xdg_from->xbd.DI1;
  xdg_to->xbd.DI1_U1 = xdg_from->xbd.DI1_U1;
  xdg_to->xbd.DI1_T1 = xdg_from->xbd.DI1_T1;
  xdg_to->xbd.DI1_D1 = xdg_from->xbd.DI1_D1;
  xdg_to->xbd.DI1_S1 = xdg_from->xbd.DI1_S1;
  xdg_to->xbd.DI1_MS = xdg_from->xbd.DI1_MS;
  xdg_to->xbd.DI1_RE = xdg_from->xbd.DI1_RE;
  xdg_to->xbd.US1 = xdg_from->xbd.US1;
  xdg_to->xbd.US1_U1 = xdg_from->xbd.US1_U1;
  xdg_to->xbd.US1_T1 = xdg_from->xbd.US1_T1;
  xdg_to->xbd.US1_D1 = xdg_from->xbd.US1_D1;
  xdg_to->xbd.US1_MS = xdg_from->xbd.US1_MS;
  xdg_to->xbd.US1_RE = xdg_from->xbd.US1_RE;
  xdg_to->xbd.CQ1 = xdg_from->xbd.CQ1;
  xdg_to->xbd.CQ1_U1 = xdg_from->xbd.CQ1_U1;
  xdg_to->xbd.CQ1_T1 = xdg_from->xbd.CQ1_T1;
  xdg_to->xbd.CQ1_D1 = xdg_from->xbd.CQ1_D1;
  xdg_to->xbd.CQ1_MS = xdg_from->xbd.CQ1_MS;
  xdg_to->xbd.CQ1_RE = xdg_from->xbd.CQ1_RE;
  xdg_to->xbd.DE1 = xdg_from->xbd.DE1;
  xdg_to->xbd.DE1_U1 = xdg_from->xbd.DE1_U1;
  xdg_to->xbd.DE1_T1 = xdg_from->xbd.DE1_T1;
  xdg_to->xbd.DE1_D1 = xdg_from->xbd.DE1_D1;
  xdg_to->xbd.DE1_MS = xdg_from->xbd.DE1_MS;
  xdg_to->xbd.X2 = xdg_from->xbd.X2;
  xdg_to->xbd.U2 = xdg_from->xbd.U2;
  xdg_to->xbd.T2 = xdg_from->xbd.T2;
  xdg_to->xbd.D2 = xdg_from->xbd.D2;
  xdg_to->xbd.S2 = xdg_from->xbd.S2;
  xdg_to->xbd.AMPL2 = xdg_from->xbd.AMPL2;
  xdg_to->xbd.U2_UEI = xdg_from->xbd.U2_UEI;
  xdg_to->xbd.U2_MS = xdg_from->xbd.U2_MS;
  xdg_to->xbd.DW2 = xdg_from->xbd.DW2;
  xdg_to->xbd.H2 = xdg_from->xbd.H2;
  xdg_to->xbd.H2_T2 = xdg_from->xbd.H2_T2;
  xdg_to->xbd.H2_D2 = xdg_from->xbd.H2_D2;
  xdg_to->xbd.M2 = xdg_from->xbd.M2;
  xdg_to->xbd.M2_U2 = xdg_from->xbd.M2_U2;
  xdg_to->xbd.M2_MS = xdg_from->xbd.M2_MS;
  xdg_to->xbd.R2 = xdg_from->xbd.R2;
  xdg_to->xbd.R2_U2 = xdg_from->xbd.R2_U2;
  xdg_to->xbd.R2_MS = xdg_from->xbd.R2_MS;
  xdg_to->xbd.V2 = xdg_from->xbd.V2;
  xdg_to->xbd.V2_U2 = xdg_from->xbd.V2_U2;
  xdg_to->xbd.V2_MS = xdg_from->xbd.V2_MS;
  xdg_to->xbd.V2_RE = xdg_from->xbd.V2_RE;
  xdg_to->xbd.HK2 = xdg_from->xbd.HK2;
  xdg_to->xbd.HK2_U2 = xdg_from->xbd.HK2_U2;
  xdg_to->xbd.HK2_T2 = xdg_from->xbd.HK2_T2;
  xdg_to->xbd.HK2_D2 = xdg_from->xbd.HK2_D2;
  xdg_to->xbd.HK2_MS = xdg_from->xbd.HK2_MS;
  xdg_to->xbd.HS2 = xdg_from->xbd.HS2;
  xdg_to->xbd.HS2_U2 = xdg_from->xbd.HS2_U2;
  xdg_to->xbd.HS2_T2 = xdg_from->xbd.HS2_T2;
  xdg_to->xbd.HS2_D2 = xdg_from->xbd.HS2_D2;
  xdg_to->xbd.HS2_MS = xdg_from->xbd.HS2_MS;
  xdg_to->xbd.HS2_RE = xdg_from->xbd.HS2_RE;
  xdg_to->xbd.HC2 = xdg_from->xbd.HC2;
  xdg_to->xbd.HC2_U2 = xdg_from->xbd.HC2_U2;
  xdg_to->xbd.HC2_T2 = xdg_from->xbd.HC2_T2;
  xdg_to->xbd.HC2_D2 = xdg_from->xbd.HC2_D2;
  xdg_to->xbd.HC2_MS = xdg_from->xbd.HC2_MS;
  xdg_to->xbd.RT2 = xdg_from->xbd.RT2;
  xdg_to->xbd.RT2_U2 = xdg_from->xbd.RT2_U2;
  xdg_to->xbd.RT2_T2 = xdg_from->xbd.RT2_T2;
  xdg_to->xbd.RT2_MS = xdg_from->xbd.RT2_MS;
  xdg_to->xbd.RT2_RE = xdg_from->xbd.RT2_RE;
  xdg_to->xbd.CF2 = xdg_from->xbd.CF2;
  xdg_to->xbd.CF2_U2 = xdg_from->xbd.CF2_U2;
  xdg_to->xbd.CF2_T2 = xdg_from->xbd.CF2_T2;
  xdg_to->xbd.CF2_D2 = xdg_from->xbd.CF2_D2;
  xdg_to->xbd.CF2_MS = xdg_from->xbd.CF2_MS;
  xdg_to->xbd.CF2_RE = xdg_from->xbd.CF2_RE;
  xdg_to->xbd.DI2 = xdg_from->xbd.DI2;
  xdg_to->xbd.DI2_U2 = xdg_from->xbd.DI2_U2;
  xdg_to->xbd.DI2_T2 = xdg_from->xbd.DI2_T2;
  xdg_to->xbd.DI2_D2 = xdg_from->xbd.DI2_D2;
  xdg_to->xbd.DI2_S2 = xdg_from->xbd.DI2_S2;
  xdg_to->xbd.DI2_MS = xdg_from->xbd.DI2_MS;
  xdg_to->xbd.DI2_RE = xdg_from->xbd.DI2_RE;
  xdg_to->xbd.US2 = xdg_from->xbd.US2;
  xdg_to->xbd.US2_U2 = xdg_from->xbd.US2_U2;
  xdg_to->xbd.US2_T2 = xdg_from->xbd.US2_T2;
  xdg_to->xbd.US2_D2 = xdg_from->xbd.US2_D2;
  xdg_to->xbd.US2_MS = xdg_from->xbd.US2_MS;
  xdg_to->xbd.US2_RE = xdg_from->xbd.US2_RE;
  xdg_to->xbd.CQ2 = xdg_from->xbd.CQ2;
  xdg_to->xbd.CQ2_U2 = xdg_from->xbd.CQ2_U2;
  xdg_to->xbd.CQ2_T2 = xdg_from->xbd.CQ2_T2;
  xdg_to->xbd.CQ2_D2 = xdg_from->xbd.CQ2_D2;
  xdg_to->xbd.CQ2_MS = xdg_from->xbd.CQ2_MS;
  xdg_to->xbd.CQ2_RE = xdg_from->xbd.CQ2_RE;
  xdg_to->xbd.DE2 = xdg_from->xbd.DE2;
  xdg_to->xbd.DE2_U2 = xdg_from->xbd.DE2_U2;
  xdg_to->xbd.DE2_T2 = xdg_from->xbd.DE2_T2;
  xdg_to->xbd.DE2_D2 = xdg_from->xbd.DE2_D2;
  xdg_to->xbd.DE2_MS = xdg_from->xbd.DE2_MS;
  xdg_to->xbd.CFM = xdg_from->xbd.CFM;
  xdg_to->xbd.CFM_MS = xdg_from->xbd.CFM_MS;
  xdg_to->xbd.CFM_RE = xdg_from->xbd.CFM_RE;
  xdg_to->xbd.CFM_U1 = xdg_from->xbd.CFM_U1;
  xdg_to->xbd.CFM_T1 = xdg_from->xbd.CFM_T1;
  xdg_to->xbd.CFM_D1 = xdg_from->xbd.CFM_D1;
  xdg_to->xbd.CFM_U2 = xdg_from->xbd.CFM_U2;
  xdg_to->xbd.CFM_T2 = xdg_from->xbd.CFM_T2;
  xdg_to->xbd.CFM_D2 = xdg_from->xbd.CFM_D2;
  xdg_to->xbd.XT = xdg_from->xbd.XT;
  xdg_to->xbd.XT_A1 = xdg_from->xbd.XT_A1;
  xdg_to->xbd.XT_MS = xdg_from->xbd.XT_MS;
  xdg_to->xbd.XT_RE = xdg_from->xbd.XT_RE;
  xdg_to->xbd.XT_XF = xdg_from->xbd.XT_XF;
  xdg_to->xbd.XT_X1 = xdg_from->xbd.XT_X1;
  xdg_to->xbd.XT_T1 = xdg_from->xbd.XT_T1;
  xdg_to->xbd.XT_D1 = xdg_from->xbd.XT_D1;
  xdg_to->xbd.XT_U1 = xdg_from->xbd.XT_U1;
  xdg_to->xbd.XT_X2 = xdg_from->xbd.XT_X2;
  xdg_to->xbd.XT_T2 = xdg_from->xbd.XT_T2;
  xdg_to->xbd.XT_D2 = xdg_from->xbd.XT_D2;
  xdg_to->xbd.XT_U2 = xdg_from->xbd.XT_U2;
  xdg_to->xbd.DWTE = xdg_from->xbd.DWTE;
  xdg_to->xbd.QINFBL = xdg_from->xbd.QINFBL;
  xdg_to->xbd.TKBL = xdg_from->xbd.TKBL;
  xdg_to->xbd.TKBL_MS = xdg_from->xbd.TKBL_MS;
  xdg_to->xbd.RSTBL = xdg_from->xbd.RSTBL;
  xdg_to->xbd.RSTBL_MS = xdg_from->xbd.RSTBL_MS;
  xdg_to->xbd.HSTINV = xdg_from->xbd.HSTINV;
  xdg_to->xbd.HSTINV_MS = xdg_from->xbd.HSTINV_MS;
  xdg_to->xbd.REYBL = xdg_from->xbd.REYBL;
  xdg_to->xbd.REYBL_MS = xdg_from->xbd.REYBL_MS;
  xdg_to->xbd.REYBL_RE = xdg_from->xbd.REYBL_RE;
  xdg_to->xbd.GAMBL = xdg_from->xbd.GAMBL;
  xdg_to->xbd.GM1BL = xdg_from->xbd.GM1BL;
  xdg_to->xbd.HVRA = xdg_from->xbd.HVRA;
  xdg_to->xbd.BULE = xdg_from->xbd.BULE;
  xdg_to->xbd.XIFORC = xdg_from->xbd.XIFORC;
  xdg_to->xbd.AMCRIT = xdg_from->xbd.AMCRIT;
  xdg_to->bld.SCCON = xdg_from->bld.SCCON;
  xdg_to->bld.GACON = xdg_from->bld.GACON;
  xdg_to->bld.GBCON = xdg_from->bld.GBCON;
  xdg_to->bld.GCCON = xdg_from->bld.GCCON;
  xdg_to->bld.DLCON = xdg_from->bld.DLCON;
  xdg_to->bld.CTRCON = xdg_from->bld.CTRCON;
  xdg_to->bld.CTRCEX = xdg_from->bld.CTRCEX;
  xdg_to->bld.DUXCON = xdg_from->bld.DUXCON;
  xdg_to->bld.CTCON = xdg_from->bld.CTCON;
  xdg_to->bld.CFFAC = xdg_from->bld.CFFAC;
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.GAM[i] = xdg_from->xfd.GAM[i];
  }
  for ( i = 0; i < IQX*2; i++ )
  {
    xdg_to->xfd.GAMU[i] = xdg_from->xfd.GAMU[i];
  }
  for ( i = 0; i < IZX*2; i++ )
  {
    xdg_to->xfd.QINVU[i] = xdg_from->xfd.QINVU[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.GAM_A[i] = xdg_from->xfd.GAM_A[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.QINV[i] = xdg_from->xfd.QINV[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.QINV_A[i] = xdg_from->xfd.QINV_A[i];
  }
  for ( i = 0; i < IQX*IQX; i++ )
  {
    xdg_to->xfd.AIJ[i] = xdg_from->xfd.AIJ[i];
  }
  for ( i = 0; i < IQX*IZX; i++ )
  {
    xdg_to->xfd.BIJ[i] = xdg_from->xfd.BIJ[i];
  }
  for ( i = 0; i < IZX*IZX; i++ )
  {
    xdg_to->xfd.DIJ[i] = xdg_from->xfd.DIJ[i];
  }
  for ( i = 0; i < IWX*IQX; i++ )
  {
    xdg_to->xfd.CIJ[i] = xdg_from->xfd.CIJ[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.DZDG[i] = xdg_from->xfd.DZDG[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.DZDN[i] = xdg_from->xfd.DZDN[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.DQDG[i] = xdg_from->xfd.DQDG[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.DZDM[i] = xdg_from->xfd.DZDM[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.DQDM[i] = xdg_from->xfd.DQDM[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.X[i] = xdg_from->xfd.X[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.Y[i] = xdg_from->xfd.Y[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.NX[i] = xdg_from->xfd.NX[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.NY[i] = xdg_from->xfd.NY[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.S[i] = xdg_from->xfd.S[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.APANEL[i] = xdg_from->xfd.APANEL[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.SIG[i] = xdg_from->xfd.SIG[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.XP[i] = xdg_from->xfd.XP[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.YP[i] = xdg_from->xfd.YP[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.QF0[i] = xdg_from->xfd.QF0[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.QF1[i] = xdg_from->xfd.QF1[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.QF2[i] = xdg_from->xfd.QF2[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.QF3[i] = xdg_from->xfd.QF3[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.AIJPIV[i] = xdg_from->xfd.AIJPIV[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.IBLTE[i] = xdg_from->xfd.IBLTE[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.NBL[i] = xdg_from->xfd.NBL[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.IPAN[i] = xdg_from->xfd.IPAN[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.ISYS[i] = xdg_from->xfd.ISYS[i];
  }
  for ( i = 0; i < IBX; i++ )
  {
    xdg_to->xfd.XB[i] = xdg_from->xfd.XB[i];
  }
  for ( i = 0; i < IBX; i++ )
  {
    xdg_to->xfd.YB[i] = xdg_from->xfd.YB[i];
  }
  for ( i = 0; i < IBX; i++ )
  {
    xdg_to->xfd.SB[i] = xdg_from->xfd.SB[i];
  }
  for ( i = 0; i < IBX; i++ )
  {
    xdg_to->xfd.XBP[i] = xdg_from->xfd.XBP[i];
  }
  for ( i = 0; i < IBX; i++ )
  {
    xdg_to->xfd.YBP[i] = xdg_from->xfd.YBP[i];
  }
  for ( i = 0; i < 5*IBX; i++ )
  {
    xdg_to->xfd.SNEW[i] = xdg_from->xfd.SNEW[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W1[i] = xdg_from->xfd.W1[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W2[i] = xdg_from->xfd.W2[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W3[i] = xdg_from->xfd.W3[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W4[i] = xdg_from->xfd.W4[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W5[i] = xdg_from->xfd.W5[i];
  }
  for ( i = 0; i < 6*IQX; i++ )
  {
    xdg_to->xfd.W6[i] = xdg_from->xfd.W6[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.CPI[i] = xdg_from->xfd.CPI[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.CPV[i] = xdg_from->xfd.CPV[i];
  }
  for ( i = 0; i < IZX; i++ )
  {
    xdg_to->xfd.QVIS[i] = xdg_from->xfd.QVIS[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.VTI[i] = xdg_from->xfd.VTI[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.XSSI[i] = xdg_from->xfd.XSSI[i];
  }
  for ( i = 0; i < IWX; i++ )
  {
    xdg_to->xfd.WGAP[i] = xdg_from->xfd.WGAP[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.XSTRIP[i] = xdg_from->xfd.XSTRIP[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.XSSITR[i] = xdg_from->xfd.XSSITR[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.UINV[i] = xdg_from->xfd.UINV[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.UINV_A[i] = xdg_from->xfd.UINV_A[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.UEDG[i] = xdg_from->xfd.UEDG[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.THET[i] = xdg_from->xfd.THET[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.DSTR[i] = xdg_from->xfd.DSTR[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.CTAU[i] = xdg_from->xfd.CTAU[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.MASS[i] = xdg_from->xfd.MASS[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.TAU[i] = xdg_from->xfd.TAU[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.DIS[i] = xdg_from->xfd.DIS[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.CTQ[i] = xdg_from->xfd.CTQ[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.DELT[i] = xdg_from->xfd.DELT[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.TSTR[i] = xdg_from->xfd.TSTR[i];
  }
  for ( i = 0; i < IVX*ISX; i++ )
  {
    xdg_to->xfd.USLP[i] = xdg_from->xfd.USLP[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.ITRAN[i] = xdg_from->xfd.ITRAN[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.TFORCE[i] = xdg_from->xfd.TFORCE[i];
  }
  for ( i = 0; i < 3*IZX*IZX; i++ )
  {
    xdg_to->xfd.VM[i] = xdg_from->xfd.VM[i];
  }
  for ( i = 0; i < 3*2*IZX; i++ )
  {
    xdg_to->xfd.VA[i] = xdg_from->xfd.VA[i];
  }
  for ( i = 0; i < 3*2*IZX; i++ )
  {
    xdg_to->xfd.VB[i] = xdg_from->xfd.VB[i];
  }
  for ( i = 0; i < 3*2*IZX; i++ )
  {
    xdg_to->xfd.VDEL[i] = xdg_from->xfd.VDEL[i];
  }
  for ( i = 0; i < 3*2; i++ )
  {
    xdg_to->xfd.VZ[i] = xdg_from->xfd.VZ[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.XOCTR[i] = xdg_from->xfd.XOCTR[i];
  }
  for ( i = 0; i < ISX; i++ )
  {
    xdg_to->xfd.YOCTR[i] = xdg_from->xfd.YOCTR[i];
  }
  for ( i = 0; i < IVX*2; i++ )
  {
    xdg_to->xfd.UNEW[i] = xdg_from->xfd.UNEW[i];
  }
  for ( i = 0; i < IVX*2; i++ )
  {
    xdg_to->xfd.U_AC[i] = xdg_from->xfd.U_AC[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.QNEW[i] = xdg_from->xfd.QNEW[i];
  }
  for ( i = 0; i < IQX; i++ )
  {
    xdg_to->xfd.Q_AC[i] = xdg_from->xfd.Q_AC[i];
  }
  for ( i = 0; i < NCOM; i++ )
  {
    xdg_to->xbd.C1SAV[i] = xdg_from->xbd.C1SAV[i];
  }
  for ( i = 0; i < NCOM; i++ )
  {
    xdg_to->xbd.C2SAV[i] = xdg_from->xbd.C2SAV[i];
  }
  for ( i = 0; i < 4*5; i++ )
  {
    xdg_to->xbd.VS1[i] = xdg_from->xbd.VS1[i];
  }
  for ( i = 0; i < 4*5; i++ )
  {
    xdg_to->xbd.VS2[i] = xdg_from->xbd.VS2[i];
  }
  for ( i = 0; i < 4; i++ )
  {
    xdg_to->xbd.VSREZ[i] = xdg_from->xbd.VSREZ[i];
  }
  for ( i = 0; i < 4; i++ )
  {
    xdg_to->xbd.VSR[i] = xdg_from->xbd.VSR[i];
  }
  for ( i = 0; i < 4; i++ )
  {
    xdg_to->xbd.VSM[i] = xdg_from->xbd.VSM[i];
  }
  for ( i = 0; i < 4; i++ )
  {
    xdg_to->xbd.VSX[i] = xdg_from->xbd.VSX[i];
  }
}
