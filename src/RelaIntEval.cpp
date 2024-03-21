#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "RelaIntEval.h"
#include <mpi.h>
using namespace std;

void RelaIntegration::MetropolisData(fstream& infoFile, fstream& basisFile, fstream& xyzFile, fstream& coeffsFile, fstream& orbEneFile) {

  rB.readBasisRel(infoFile, basisFile, xyzFile, coeffsFile, orbEneFile);
  lambda = 2.0 * (rB.bI.orbEne[rB.bI.occOrb] - rB.bI.orbEne[rB.bI.occOrb-1]);
	count = MRC.numberWalkers*(MRC.numberWalkers - 1) * 0.5;
  
   for (int i=0; i<2; i++) {
    occs[i].factor.resize(256);
    virsA[i].factor.resize(256);
    virsB[i].factor.resize(256);
    vir13[i].factor.resize(16);
    vir14[i].factor.resize(16);
    vir23[i].factor.resize(16);
    vir24[i].factor.resize(16);
    occ13[i].factor.resize(16);
    occ24[i].factor.resize(16);
  }

}

void RelaIntegration::RelaMP2Energy (int mynode, int totalnodes, fstream& infoFile, fstream& basisFile, fstream& xyzFile, fstream& coeffsFile, fstream& orbEneFile, fstream& MCFile) {

	double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, r12, r34, w12, w34;
	double emp2a, emp2b, MP2Ener, scsEMP2, result, denom;
	double diff, emp2Ave, scsEMP2Ave, emp2Var, scsEMP2Var, Rvalue, wgTau, tau;
	unsigned int tot;
  double MCintegralA, MCintegralB, movLen, Eg;
  int h, startval, endval, nS, nW, local_offset, rows_local, Nb, Nbs, ni;
  MPI_Status status;

	MRC.Data(xyzFile, MCFile);
  rB.readBasisRel(infoFile, basisFile, xyzFile, coeffsFile, orbEneFile);
  lambda = 2.0 * (rB.bI.orbEne[rB.bI.occOrb] - rB.bI.orbEne[rB.bI.occOrb-1]);
  count = MRC.numberWalkers*(MRC.numberWalkers - 1) * 0.5;

   for (int i=0; i<2; i++) {
    occs[i].factor.resize(256);
    virsA[i].factor.resize(256);
    virsB[i].factor.resize(256);
    vir13[i].factor.resize(16);
    vir14[i].factor.resize(16);
    vir23[i].factor.resize(16);
    vir24[i].factor.resize(16);
    occ13[i].factor.resize(16);
    occ24[i].factor.resize(16);
  }

/*  
	x1 = 0.1;
	y1 = 0.2;
	z1 = 0.3;
	cout << "CHECKING: " << endl;
	rB.mofunc(x1, y1, z1, 0);
*/

	tot = 0;
	emp2Ave = 0.0;
	emp2Var = 0.0;
	FinalMP2Ener = 0.0;
  nS = MRC.nSteps;
  nW = MRC.numberWalkers;
  Eg = MRC.Eg;
  Nbs = MRC.Nbsize;
  Nb = nS/Nbs;

  ////////////////STARTS MPI


  if (mynode == 0) {
    for (int j=1; j<totalnodes; j++) {
      MPI_Send(&nS,1,MPI_INT,j,0,MPI_COMM_WORLD);
      MPI_Send(&nW,1,MPI_INT,j,0,MPI_COMM_WORLD);
      MPI_Send(&Eg,1,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
     }
  } else
    if (mynode != 0) { 
      MPI_Recv(&nS,1,MPI_INT,0,0,MPI_COMM_WORLD, &status);
      MPI_Recv(&nW,1,MPI_INT,0,0,MPI_COMM_WORLD, &status);
      MPI_Recv(&Eg,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD, &status);
    }

  startval = nS*mynode/totalnodes;
  endval = nS*(mynode + 1)/totalnodes - 1;

  rows_local = (int) floor(nS/totalnodes);
  local_offset = mynode*rows_local;
  if(mynode == (totalnodes-1))
    rows_local = nS - rows_local*(totalnodes-1);

  MRC.motion(mynode);
 
  walker.random();
  ni = 0;
  for (int ij=0; ij<Nb; ij++) {
   if (Nb != nS) {
     tot = 0;
     emp2Ave = 0.0;
     emp2Var = 0.0;
     FinalMP2Ener = 0.0;
   }
	  for (int n=0; n<Nbs; n++) {
		  MP2Ener = 0.0;
		  Rvalue = walker.uniform();
		  wgTau = 1.0 / (lambda * (1.0 - Rvalue));
		  tau = -log(1.0-Rvalue)/lambda;
		  for (int i=0; i<MRC.numberWalkers-1; i++) {
			
			  x1 = MRC.wp[i].rSteps[ni].e1[0];
        y1 = MRC.wp[i].rSteps[ni].e1[1];
        z1 = MRC.wp[i].rSteps[ni].e1[2];
			  rB.mofunc(x1, y1, z1, 0);

        x2 = MRC.wp[i].rSteps[ni].e2[0];
        y2 = MRC.wp[i].rSteps[ni].e2[1];
        z2 = MRC.wp[i].rSteps[ni].e2[2];
			  rB.mofunc(x2, y2, z2, 1);			

        r12 = MRC.wp[i].rSteps[ni].r12;
        w12 = MRC.wp[i].rSteps[ni].w12;

			  for (int j=i+1; j<MRC.numberWalkers; j++) {

          x3 = MRC.wp[j].rSteps[ni].e1[0];
          y3 = MRC.wp[j].rSteps[ni].e1[1];
          z3 = MRC.wp[j].rSteps[ni].e1[2];
				  rB.mofunc(x3, y3, z3, 2);			

          x4 = MRC.wp[j].rSteps[ni].e2[0];
          y4 = MRC.wp[j].rSteps[ni].e2[1];
          z4 = MRC.wp[j].rSteps[ni].e2[2];
				  rB.mofunc(x4, y4, z4, 3);			

          r34 = MRC.wp[j].rSteps[ni].r12;
       	  w34 = MRC.wp[j].rSteps[ni].w12;

			    denom = 1.0/(r12*r34*w12*w34);

			    oIJTau(occ13, tau, 0, 2);	
			    oIJTau(occ24, tau, 1, 3);
			    vABTau(vir13, tau, 0, 2);
			    vABTau(vir24, tau, 1, 3);
			    vABTau(vir23, tau, 1, 2);
			    vABTau(vir14, tau, 0, 3);

          h = 0;
          for (int m=0; m<16; m++) {
            for (int n=0; n<16; n++) {
              occs[0].factor[h] = rB.complexReal(occ13[0].factor[m], occ13[1].factor[m], occ24[0].factor[n], occ24[1].factor[n]);
              occs[1].factor[h] = rB.complexIm(occ13[0].factor[m], occ13[1].factor[m], occ24[0].factor[n], occ24[1].factor[n]);
              virsA[0].factor[h] = rB.complexReal(vir13[0].factor[m], vir13[1].factor[m], vir24[0].factor[n], vir24[1].factor[n]);
              virsA[1].factor[h] = rB.complexIm(vir13[0].factor[m], vir13[1].factor[m], vir24[0].factor[n], vir24[1].factor[n]);
              virsB[0].factor[h] = rB.complexReal(vir14[0].factor[m], vir14[1].factor[m], vir23[0].factor[n], vir23[1].factor[n]);
              virsB[1].factor[h] = rB.complexIm(vir23[0].factor[m], vir23[1].factor[m], vir14[0].factor[n], vir14[1].factor[n]);
              h++;
            }
          }

          MCintegralA = 0.0;
          MCintegralB = 0.0;
          for (int k=0; k<256; k++) {
            MCintegralA += rB.complexReal(occs[0].factor[k], occs[1].factor[k], virsA[0].factor[k], virsA[1].factor[k]);
            MCintegralB += rB.complexReal(occs[0].factor[k], occs[1].factor[k], virsB[0].factor[k], virsB[1].factor[k]);
          }

          MCintegralA *= wgTau;
          MCintegralB *= wgTau;
          emp2a = -0.5*abs(MCintegralA-MCintegralB)*denom;
	
			    MP2Ener += emp2a;
			  }

		  }

		  result = MP2Ener/count;

		  diff = result - emp2Ave;

      if (abs(diff) <= 100.0) {
		    tot++;				
		    FinalMP2Ener += result;
  		  emp2Ave += diff/tot;
  		  emp2Var += (diff*diff*double((tot-1.0)/tot));
		  //  cout << setprecision(36) << ni << "   " << result << "\t " << diff*diff*double((tot-1.0)/tot) << endl;
        }
      ni++;
      }
		  cout << setprecision(36) << FinalMP2Ener/tot << "\t " << emp2Var << endl;
    }

	FinalMP2Ener /= tot;

//  result = sqrt(emp2Var/double(tot*(tot-1.0)));

//	cout << "The result in the node: " << mynode << " is  " << FinalMP2Ener <<  " with final uncertainty:    " << result << endl;

}

void RelaIntegration::oIJ(double &aval, double &bval, int t, int ik, int jk) {

	aval = 0.0;
	bval = 0.0;

	for (int i=0; i<rB.bI.occOrb; i++) 
		for (int m=0; m<4; m++) {
			aval += wg * rB.complexReal(rB.bI.moSpi[ik].Occvalmo[m].valuesMos[0][i], rB.bI.moSpi[ik].Occvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[0][i], -rB.bI.moSpi[jk].Occvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].occExpo[i];
			bval += wg * rB.complexIm(rB.bI.moSpi[ik].Occvalmo[m].valuesMos[0][i], rB.bI.moSpi[ik].Occvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[0][i], -rB.bI.moSpi[jk].Occvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].occExpo[i];
		}

}

void RelaIntegration::vABTau(array<values, 2>& vir, double tau, int ik, int jk) {

  double tauValue;
  int h = 0;
  for (int j=0; j<4; j++)
    for (int k=0; k<4; k++) {
      vir[0].factor[h] = 0.0;
      vir[1].factor[h] = 0.0;
      h++;
    }

  for (int i=0; i<rB.bI.virOrb; i++) {
    tauValue = exp(-tau * rB.bI.orbEne[rB.bI.occOrb + i]);
    h = 0;
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++) {
        vir[0].factor[h] += wg * rB.complexReal(rB.bI.moSpi[ik].Virvalmo[j].valuesMos[0][i], rB.bI.moSpi[ik].Virvalmo[j].valuesMos[1][i], rB.bI.moSpi[jk].Virvalmo[k].valuesMos[0][i], -rB.bI.moSpi[jk].Virvalmo[k].valuesMos[1][i]) * tauValue;;
        vir[1].factor[h] += wg * rB.complexIm(rB.bI.moSpi[ik].Virvalmo[j].valuesMos[0][i], rB.bI.moSpi[ik].Virvalmo[j].valuesMos[1][i], rB.bI.moSpi[jk].Virvalmo[k].valuesMos[0][i], -rB.bI.moSpi[jk].Virvalmo[k].valuesMos[1][i]) * tauValue;
        h++;
      }
  }

}

void RelaIntegration::oIJTau(array<values, 2>&  occ, double tau, int ik, int jk) {

  double tauValue;
  int h = 0;
  for (int j=0; j<4; j++)
    for (int k=0; k<4; k++) {
      occ[0].factor[h] = 0.0;
      occ[1].factor[h] = 0.0;
      h++;
    }


  for (int i=0; i<rB.bI.occOrb; i++) {
		tauValue = exp(tau * rB.bI.orbEne[i]);
    h = 0;
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++) {
        occ[0].factor[h] += wg * rB.complexReal(rB.bI.moSpi[ik].Occvalmo[j].valuesMos[0][i], rB.bI.moSpi[ik].Occvalmo[j].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[k].valuesMos[0][i], -rB.bI.moSpi[jk].Occvalmo[k].valuesMos[1][i]) * tauValue;
        occ[1].factor[h] += wg * rB.complexIm(rB.bI.moSpi[ik].Occvalmo[j].valuesMos[0][i], rB.bI.moSpi[ik].Occvalmo[j].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[k].valuesMos[0][i], -rB.bI.moSpi[jk].Occvalmo[k].valuesMos[1][i]) * tauValue;
        h++;
      }
  }
}

void RelaIntegration::oIJK(double &aval, double &bval, int t, int ik, int jk) {

	aval = 0.0;
	bval = 0.0;

	for (int i=0; i<rB.bI.occOrb; i++) 
		for (int m=0; m<4; m++) {
			aval += wg * rB.complexReal(rB.bI.moSpi[ik].Occvalmo[m].valuesMos[0][i], -rB.bI.moSpi[ik].Occvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[0][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].occExpo[i];
			bval += wg * rB.complexIm(rB.bI.moSpi[ik].Occvalmo[m].valuesMos[0][i], -rB.bI.moSpi[ik].Occvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[0][i], rB.bI.moSpi[jk].Occvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].occExpo[i];
		}
}

void RelaIntegration::vABK(double &aval, double &bval, int t, int ik, int jk) {

	aval = 0.0;
	bval = 0.0;

	for (int i=0; i<rB.bI.virOrb; i++) 
		for (int m=0; m<4; m++) {
			aval += wg * rB.complexReal(rB.bI.moSpi[ik].Virvalmo[m].valuesMos[0][i], rB.bI.moSpi[ik].Virvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Virvalmo[m].valuesMos[0][i], -rB.bI.moSpi[jk].Virvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].virExpo[i];
			bval += wg * rB.complexIm(rB.bI.moSpi[ik].Virvalmo[m].valuesMos[0][i], rB.bI.moSpi[ik].Virvalmo[m].valuesMos[1][i], rB.bI.moSpi[jk].Virvalmo[m].valuesMos[0][i], -rB.bI.moSpi[jk].Virvalmo[m].valuesMos[1][i]) * RelaTauExpos[t].virExpo[i];
		}

}

void RelaIntegration::gaussKronrodRela () {

  double tau1, tau2, fac;
	int v1;

  for (int i=0; i<10; i++) {
    tau1 = 0.5 - 0.5*gk.tj[i];
    tau2 = 0.5 + 0.5*gk.tj[i];
    gk.tau[2*i] = tau1;
    gk.tau[2*i+1] = tau2;
    gk.tp[2*i] = (1.0 - tau1)/tau1;
    gk.tp[2*i+1] = (1.0 - tau2)/tau2;
  }
    gk.tau[20] = 0.5;
    gk.tp[20] = 1.0;
	
	v1 = rB.bI.occOrb;
  for (int t=0; t<21; t++) {
    RelaTauExpos[t].occExpo.resize(rB.bI.occOrb);
		RelaTauExpos[t].virExpo.resize(rB.bI.virOrb);
    for (int i=0; i<rB.bI.occOrb; i++) {  
      RelaTauExpos[t].occExpo[i] = exp(gk.tp[t]*rB.bI.orbEne[i]);
	//		cout << i << "   " << gk.tp[t]*rB.bI.orbEne[i] << "  =  " << RelaTauExpos[t].occExpo[i] << endl;
		}
		for (int j=0; j<rB.bI.virOrb; j++) {
      RelaTauExpos[t].virExpo[j] = exp(-gk.tp[t]*rB.bI.orbEne[v1+j]);
	//		cout << j << "   " << gk.tp[t]*rB.bI.orbEne[v1+j] << "  =  " << RelaTauExpos[t].virExpo[j] << endl;
		}

    gk.denom[t] = 0.5*gk.wj[t]/(gk.tau[t]*gk.tau[t]);
  }
}

void RelaIntegration::tauIntegration() {

	double Rvalue;
	int k1;
	walker.random();
  lambda = 2.0 * (rB.bI.orbEne[rB.bI.occOrb] - rB.bI.orbEne[rB.bI.occOrb-1]);
	wgTaus.resize(MRC.nSteps);
	taus.resize(MRC.nSteps);

	for (int i=0; i<MRC.nSteps; i++) {
		Rvalue = walker.uniform();
		wgTaus[i] = 1.0 / (lambda * (1.0 - Rvalue));
		taus[i] = -log(1.0-Rvalue)/lambda;
	}


}

