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

#include "dtype.h"

//int nclv;      // number of microphysics variables
//int ncldql;    // liquid cloud water
//int ncldqi;    // ice cloud water
//int ncldqr;    // rain water
//int ncldqs;    // snow
//int ncldqv;    // vapour

struct TECLDP {
  dtype ramid;
  dtype rcldiff;
  dtype rcldiff_convi;
  dtype rclcrit;
  dtype rclcrit_sea;
  dtype rclcrit_land;
  dtype rkconv;
  dtype rprc1;
  dtype rprc2;
  dtype rcldmax;
  dtype rpecons;
  dtype rvrfactor;
  dtype rprecrhmax;
  dtype rtaumel;
  dtype ramin;
  dtype rlmin;
  dtype rkooptau;
  dtype rcldtopp;
  dtype rlcritsnow;
  dtype rsnowlin1;
  dtype rsnowlin2;
  dtype ricehi1;
  dtype ricehi2;
  dtype riceinit;
  dtype rvice;
  dtype rvrain;
  dtype rvsnow;
  dtype rthomo;
  dtype rcovpmin;
  dtype rccn;
  dtype rnice;
  dtype rccnom;
  dtype rccnss;
  dtype rccnsu;
  dtype rcldtopcf;
  dtype rdepliqrefrate;
  dtype rdepliqrefdepth;
  dtype rcl_kkaac;
  dtype rcl_kkbac;
  dtype rcl_kkaau;
  dtype rcl_kkbauq;
  dtype rcl_kkbaun;
  dtype rcl_kk_cloud_num_sea;
  dtype rcl_kk_cloud_num_land;
  dtype rcl_ai;
  dtype rcl_bi;
  dtype rcl_ci;
  dtype rcl_di;
  dtype rcl_x1i;
  dtype rcl_x2i;
  dtype rcl_x3i;
  dtype rcl_x4i;
  dtype rcl_const1i;
  dtype rcl_const2i;
  dtype rcl_const3i;
  dtype rcl_const4i;
  dtype rcl_const5i;
  dtype rcl_const6i;
  dtype rcl_apb1;
  dtype rcl_apb2;
  dtype rcl_apb3;
  dtype rcl_as;
  dtype rcl_bs;
  dtype rcl_cs;
  dtype rcl_ds;
  dtype rcl_x1s;
  dtype rcl_x2s;
  dtype rcl_x3s;
  dtype rcl_x4s;
  dtype rcl_const1s;
  dtype rcl_const2s;
  dtype rcl_const3s;
  dtype rcl_const4s;
  dtype rcl_const5s;
  dtype rcl_const6s;
  dtype rcl_const7s;
  dtype rcl_const8s;
  dtype rdenswat;
  dtype rdensref;
  dtype rcl_ar;
  dtype rcl_br;
  dtype rcl_cr;
  dtype rcl_dr;
  dtype rcl_x1r;
  dtype rcl_x2r;
  dtype rcl_x4r;
  dtype rcl_ka273;
  dtype rcl_cdenom1;
  dtype rcl_cdenom2;
  dtype rcl_cdenom3;
  dtype rcl_schmidt;
  dtype rcl_dynvisc;
  dtype rcl_const1r;
  dtype rcl_const2r;
  dtype rcl_const3r;
  dtype rcl_const4r;
  dtype rcl_fac1;
  dtype rcl_fac2;
  dtype rcl_const5r;
  dtype rcl_const6r;
  dtype rcl_fzrab;
  dtype rcl_fzrbb;
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
  dtype nshapep;
  dtype nshapeq;
  int nbeta;
  //dtype rbeta[0][100];
  //dtype rbetap1[0][100];
} ;

//struct TECLDP *yrecldp;

#endif
