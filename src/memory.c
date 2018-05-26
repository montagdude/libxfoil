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
// Checks that xfoil_init has been called via status of xdg->xfd.GAM
//
/******************************************************************************/
void xdg_allocated(const xfoil_data_group *xdg, int *stat)
{
  if (! xdg->xfd.GAM)
    *stat = 1;
  else
    *stat = 0;
}
