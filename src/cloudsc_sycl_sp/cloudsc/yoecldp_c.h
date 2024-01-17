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

//int nclv;      // number of microphysics variables
//int ncldql;    // liquid cloud water
//int ncldqi;    // ice cloud water
//int ncldqr;    // rain water
//int ncldqs;    // snow
//int ncldqv;    // vapour

struct TECLDP {
  float ramid;
  float rcldiff;
  float rcldiff_convi;
  float rclcrit;
  float rclcrit_sea;
  float rclcrit_land;
  float rkconv;
  float rprc1;
  float rprc2;
  float rcldmax;
  float rpecons;
  float rvrfactor;
  float rprecrhmax;
  float rtaumel;
  float ramin;
  float rlmin;
  float rkooptau;
  float rcldtopp;
  float rlcritsnow;
  float rsnowlin1;
  float rsnowlin2;
  float ricehi1;
  float ricehi2;
  float riceinit;
  float rvice;
  float rvrain;
  float rvsnow;
  float rthomo;
  float rcovpmin;
  float rccn;
  float rnice;
  float rccnom;
  float rccnss;
  float rccnsu;
  float rcldtopcf;
  float rdepliqrefrate;
  float rdepliqrefdepth;
  float rcl_kkaac;
  float rcl_kkbac;
  float rcl_kkaau;
  float rcl_kkbauq;
  float rcl_kkbaun;
  float rcl_kk_cloud_num_sea;
  float rcl_kk_cloud_num_land;
  float rcl_ai;
  float rcl_bi;
  float rcl_ci;
  float rcl_di;
  float rcl_x1i;
  float rcl_x2i;
  float rcl_x3i;
  float rcl_x4i;
  float rcl_const1i;
  float rcl_const2i;
  float rcl_const3i;
  float rcl_const4i;
  float rcl_const5i;
  float rcl_const6i;
  float rcl_apb1;
  float rcl_apb2;
  float rcl_apb3;
  float rcl_as;
  float rcl_bs;
  float rcl_cs;
  float rcl_ds;
  float rcl_x1s;
  float rcl_x2s;
  float rcl_x3s;
  float rcl_x4s;
  float rcl_const1s;
  float rcl_const2s;
  float rcl_const3s;
  float rcl_const4s;
  float rcl_const5s;
  float rcl_const6s;
  float rcl_const7s;
  float rcl_const8s;
  float rdenswat;
  float rdensref;
  float rcl_ar;
  float rcl_br;
  float rcl_cr;
  float rcl_dr;
  float rcl_x1r;
  float rcl_x2r;
  float rcl_x4r;
  float rcl_ka273;
  float rcl_cdenom1;
  float rcl_cdenom2;
  float rcl_cdenom3;
  float rcl_schmidt;
  float rcl_dynvisc;
  float rcl_const1r;
  float rcl_const2r;
  float rcl_const3r;
  float rcl_const4r;
  float rcl_fac1;
  float rcl_fac2;
  float rcl_const5r;
  float rcl_const6r;
  float rcl_fzrab;
  float rcl_fzrbb;
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
  float nshapep;
  float nshapeq;
  int nbeta;
  //float rbeta[0][100];
  //float rbetap1[0][100];
} ;

//struct TECLDP *yrecldp;

#endif
