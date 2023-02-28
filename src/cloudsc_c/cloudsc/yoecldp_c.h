/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef YOECLDP_H
#define YOECLDP_H

extern int nclv;      // number of microphysics variables
extern int ncldql;    // liquid cloud water
extern int ncldqi;    // ice cloud water
extern int ncldqr;    // rain water
extern int ncldqs;    // snow
extern int ncldqv;    // vapour

struct TECLDP {
  double ramid;
  double rcldiff;
  double rcldiff_convi;
  double rclcrit;
  double rclcrit_sea;
  double rclcrit_land;
  double rkconv;
  double rprc1;
  double rprc2;
  double rcldmax;
  double rpecons;
  double rvrfactor;
  double rprecrhmax;
  double rtaumel;
  double ramin;
  double rlmin;
  double rkooptau;
  double rcldtopp;
  double rlcritsnow;
  double rsnowlin1;
  double rsnowlin2;
  double ricehi1;
  double ricehi2;
  double riceinit;
  double rvice;
  double rvrain;
  double rvsnow;
  double rthomo;
  double rcovpmin;
  double rccn;
  double rnice;
  double rccnom;
  double rccnss;
  double rccnsu;
  double rcldtopcf;
  double rdepliqrefrate;
  double rdepliqrefdepth;
  double rcl_kkaac;
  double rcl_kkbac;
  double rcl_kkaau;
  double rcl_kkbauq;
  double rcl_kkbaun;
  double rcl_kk_cloud_num_sea;
  double rcl_kk_cloud_num_land;
  double rcl_ai;
  double rcl_bi;
  double rcl_ci;
  double rcl_di;
  double rcl_x1i;
  double rcl_x2i;
  double rcl_x3i;
  double rcl_x4i;
  double rcl_const1i;
  double rcl_const2i;
  double rcl_const3i;
  double rcl_const4i;
  double rcl_const5i;
  double rcl_const6i;
  double rcl_apb1;
  double rcl_apb2;
  double rcl_apb3;
  double rcl_as;
  double rcl_bs;
  double rcl_cs;
  double rcl_ds;
  double rcl_x1s;
  double rcl_x2s;
  double rcl_x3s;
  double rcl_x4s;
  double rcl_const1s;
  double rcl_const2s;
  double rcl_const3s;
  double rcl_const4s;
  double rcl_const5s;
  double rcl_const6s;
  double rcl_const7s;
  double rcl_const8s;
  double rdenswat;
  double rdensref;
  double rcl_ar;
  double rcl_br;
  double rcl_cr;
  double rcl_dr;
  double rcl_x1r;
  double rcl_x2r;
  double rcl_x4r;
  double rcl_ka273;
  double rcl_cdenom1;
  double rcl_cdenom2;
  double rcl_cdenom3;
  double rcl_schmidt;
  double rcl_dynvisc;
  double rcl_const1r;
  double rcl_const2r;
  double rcl_const3r;
  double rcl_const4r;
  double rcl_fac1;
  double rcl_fac2;
  double rcl_const5r;
  double rcl_const6r;
  double rcl_fzrab;
  double rcl_fzrbb;
  int lcldextra, lcldbudget;
  int nssopt;
  int ncldtop;
  int naeclbc, naecldu, naeclom, naeclss, naeclsu;
  int nclddiag;
  int naercld;
  int laerliqautolsp;
  int laerliqautocp;
  int laerliqautocpb;
  int laerliqcoll;
  int laericesed;
  int laericeauto;
  double nshapep;
  double nshapeq;
  int nbeta;
  double rbeta[0][100];
  double rbetap1[0][100];
} ;

extern struct TECLDP *yrecldp;

#endif
