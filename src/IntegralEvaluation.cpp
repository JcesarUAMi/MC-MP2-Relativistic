#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include "IntegralEvaluation.h"
using namespace std;

void Integration::Data(fstream& basisFile, fstream& xyzFile, fstream& movecsFile) {

  readData.coefficientsFinal(basisFile, xyzFile, movecsFile);
  readData.nucleusBasis();
  readData.getNumOccOrb();
  MRT2.Data(xyzFile);
  MRT2.motion();
//  gaussKronrodEval();
}


void Integration::MP2Energy (fstream& xyzFile) {

	double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, r12, r34, w12, w34;
	double aval, bval, aresk, bresk, emp2a, emp2b, MP2Ener, count, scsEMP2;
	double val1, val2, result, Eg;
	double diff, emp2Ave, scsEMP2Ave, emp2Var, scsEMP2Var;
  int k;
  double Rvalue, wgTau, tau, tauValue;

	count = MRT2.numberWalkers*(MRT2.numberWalkers - 1) * 0.5;

  double lambda = 2.0 * (readData.Vir[0] - readData.Occ[readData.numberOfOccOrb-1]);

//  cout << "pos cuanto es:   " << lambda << endl;
  //  wg.Data(xyzFile);

	gfs1.resize(readData.numberOfMolOrb);
	gfs2.resize(readData.numberOfMolOrb);
	gfs3.resize(readData.numberOfMolOrb);
	gfs4.resize(readData.numberOfMolOrb);
  MosWP1.resize(readData.numberOfMolOrb); 
	MosWP2.resize(readData.numberOfMolOrb); 
	MosWP3.resize(readData.numberOfMolOrb); 
	MosWP4.resize(readData.numberOfMolOrb); 
/*
  x1 = 1.3231683558312111;
  y1 = 4.6826899164588043E-002;
  z1 = -1.2709459024724532;
	generateMoisSpherical(gfs1, x1, y1, z1);

  x2 = -2.1721055172263095;
  y2 = 1.5285385987363842;
  z2 = 0.24612123518020335;
	generateMoisSpherical(gfs2, x2, y2, z2);

  x3 = -0.17320430386932784;  
  y3 = -0.58238472897288873;
  z3 = -0.29907011179266174;
	generateMoisSpherical(gfs3, x3, y3, z3);

  x4 = 2.1982931960624068;
  y4 = 2.5675828031260783E-003;
  z4 = 0.57422224876620998;
	generateMoisSpherical(gfs4, x4, y4, z4);
  
  for (int i=0; i<readData.numberOfMolOrb; i++) {
    k = 0;
    MosWP1[i] = 0.0;
    for (int j=0; j<readData.numberOfMolOrb; j++) {
      MosWP1[i] += gfs1[j] * readData.moisCoeffs[i][j];
      MosWP2[i] += gfs2[j] * readData.moisCoeffs[i][j];
      MosWP3[i] += gfs3[j] * readData.moisCoeffs[i][j];
      MosWP4[i] += gfs4[j] * readData.moisCoeffs[i][j];
      cout << gfs1[j] << "   " << readData.moisCoeffs[i][j] << endl;
  //    cout << gfs2[j] << "   " << readData.moisCoeffs[i][j] << endl;
  //    cout << gfs3[j] << "   " << readData.moisCoeffs[i][j] << endl;
  //    cout << gfs4[j] << "   " << readData.moisCoeffs[i][j] << endl;
    }
    cout << "CAMBIO: " << i << "  " << MosWP1[i] << endl;
  //  cout << "CAMBIO: " << i << "  " << MosWP2[i] << endl;
  //  cout << "CAMBIO: " << i << "  " << MosWP3[i] << endl;
  //  cout << "CAMBIO: " << i << "  " << MosWP4[i] << endl;
  }


  aresk = 0.0;
  bresk = 0.0;

  for (int t=0; t<21; t++) {
    tauIntegral (aval, bval, t);
    aresk += gk.wj[t]*aval;
    bresk += gk.wj[t]*bval;
  }

  cout << "achi: " << aresk << "  " << bresk << endl;
  r12 = 0.0;
  r12 = (x1-x2)*(x1-x2);
  r12 += (y1-y2)*(y1-y2);
  r12 += (z1-z2)*(z1-z2);
  r12 = sqrt(r12);

  r34 = 0.0;
  r34 = (x3-x4)*(x3-x4);
  r34 += (y3-y4)*(y3-y4);
  r34 += (z3-z4)*(z3-z4);
  r34 = sqrt(r34);

  Eg = wg.constant();
  w12 = wg.weightFunction(x1, y1, z1, x2, y2, z2, r12)/Eg; 
  w34 = wg.weightFunction(x3, y3, z3, x4, y4, z4, r34)/Eg;

  cout << "ANDALEE: " << w12*w34 << "    " << r12*r34 << endl; 

  emp2a = 0.5*aresk/(r12*r34)/(w12*w34);
  emp2b = 0.5*bresk/(r12*r34)/(w12*w34);

  cout << "ANDALEE PUTITO: " << emp2a << "   " << emp2b << endl; 
	val1 = -5.0/3.0;
	val2 = 1.0/3.0;

  emp2Ave = -2.0*emp2a + emp2b;

  cout << "CHallle: " << emp2Ave << endl;
*/
	FinalMP2Ener = 0.0;
	FinalSCSEMP2 = 0.0;
	emp2Ave = 0.0;
	emp2Var = 0.0;
	scsEMP2Ave = 0.0;
	scsEMP2Var = 0.0;
	
  walker.random();  
	for (int n=0; n<MRT2.nSteps; n++) {
		MP2Ener = 0.0;
		scsEMP2 = 0.0;
    Rvalue = walker.uniform();
    wgTau = (lambda * (1.0 - Rvalue));
    tau = -log(1.0-Rvalue)/lambda;
		for (int i=0; i<MRT2.numberWalkers-1; i++) {
				
			x1 = MRT2.wp[i].rSteps[n].e1[0];
      y1 = MRT2.wp[i].rSteps[n].e1[1];
      z1 = MRT2.wp[i].rSteps[n].e1[2];
			generateMoisSpherical(gfs1, x1, y1, z1);

      x2 = MRT2.wp[i].rSteps[n].e2[0];
      y2 = MRT2.wp[i].rSteps[n].e2[1];
      z2 = MRT2.wp[i].rSteps[n].e2[2];
      r12 = MRT2.wp[i].rSteps[n].r12;
      w12 = MRT2.wp[i].rSteps[n].w12;
      generateMoisSpherical(gfs2, x2, y2, z2);
		
      for (int i=0; i<readData.numberOfMolOrb; i++) {
        k = 0;
        MosWP1[i] = 0.0;
        MosWP2[i] = 0.0;
        for (int j=0; j<readData.numberOfMolOrb; j++) {
          MosWP1[i] += gfs1[j] * readData.moisCoeffs[i][j];
          MosWP2[i] += gfs2[j] * readData.moisCoeffs[i][j];
        }
      }

			for (int j=i+1; j<MRT2.numberWalkers; j++) {

        x3 = MRT2.wp[j].rSteps[n].e1[0];
        y3 = MRT2.wp[j].rSteps[n].e1[1];
        z3 = MRT2.wp[j].rSteps[n].e1[2];
				generateMoisSpherical(gfs3, x3, y3, z3);

        x4 = MRT2.wp[j].rSteps[n].e2[0];
        y4 = MRT2.wp[j].rSteps[n].e2[1];
        z4 = MRT2.wp[j].rSteps[n].e2[2];
				generateMoisSpherical(gfs4, x4, y4, z4);

          for (int i=0; i<readData.numberOfMolOrb; i++) {
            k = 0;
            MosWP3[i] = 0.0;
            MosWP4[i] = 0.0;
            for (int j=0; j<readData.numberOfMolOrb; j++) {
              MosWP3[i] += gfs3[j] * readData.moisCoeffs[i][j];
              MosWP4[i] += gfs4[j] * readData.moisCoeffs[i][j];
            }
          }

        stochasticTauIntegration(aresk, bresk, tau);

        r34 = MRT2.wp[j].rSteps[n].r12;
       	w34 = MRT2.wp[j].rSteps[n].w12;

        /*
				for (int t=0; t<21; t++) {
					tauIntegral (aval, bval, t);	
					aresk += gk.wj[t]*aval;
					bresk += gk.wj[t]*bval;
				}
*/

				emp2a = aresk/(r12*r34)/(wgTau*w12*w34);
				emp2b = bresk/(r12*r34)/(wgTau*w12*w34);

				MP2Ener += (-2.0*emp2a + emp2b);		
				scsEMP2 += emp2a + emp2b;
		}

	}

//		result = 5.0*MP2Ener/(2.0*count);	
		result = MP2Ener/count;
			
		FinalMP2Ener += result;

    if (result >= 1000) 
		  cout << "voy n: " << n << "   " << result << "   " << FinalMP2Ener << endl;

		diff = result - emp2Ave;
		emp2Ave += diff/(n+1);
		emp2Var += (diff*diff*double(((n+1)-1.0)/(n+1)));

//		cout << "chin: " << n << "   " << diff << "  " << emp2Ave << "  " << emp2Var << endl;

	}

	FinalMP2Ener /= MRT2.nSteps;
	FinalSCSEMP2 /= MRT2.nSteps;

	result = sqrt(emp2Var/(MRT2.nSteps*(MRT2.nSteps-1.0)));

	cout << "MP2 Energy is " << FinalMP2Ener << endl;
	cout << "with variance: " << emp2Ave << "  " << result << endl;


}


void Integration::tauIntegral (double &aval, double &bval, int t) {

	int v1;
	double o13, o14, o24, o23, v13, v14, v23, v24, factor;
  o13 = 0.0;
  o14 = 0.0;
  o23 = 0.0;
  o24 = 0.0;
  v13 = 0.0;
  v14 = 0.0;
  v23 = 0.0;
  v24 = 0.0;

	factor = 0.5*gk.denom[t];
	for (int u=0; u<readData.numberOfOccOrb; u++) {
  	o13 += MosWP1[u] * MosWP3[u] * tauExpos[t].occExpo[u];
    o14 += MosWP1[u] * MosWP4[u] * tauExpos[t].occExpo[u];
    o23 += MosWP2[u] * MosWP3[u] * tauExpos[t].occExpo[u];
    o24 += MosWP2[u] * MosWP4[u] * tauExpos[t].occExpo[u];
  }

  v1 = readData.numberOfOccOrb;

  for (int v=0; v<readData.numberOfVirOrb; v++) {
    v13 += MosWP1[v1] * MosWP3[v1] * tauExpos[t].virExpo[v];
    v14 += MosWP1[v1] * MosWP4[v1] * tauExpos[t].virExpo[v];
    v23 += MosWP2[v1] * MosWP3[v1] * tauExpos[t].virExpo[v];
    v24 += MosWP2[v1] * MosWP4[v1] * tauExpos[t].virExpo[v];
    v1++;
  }
	aval = factor*(o13*o24*v13*v24 + o14*o23*v14*v23);
  bval = factor*(o14*o23*v13*v24 + o13*o24*v14*v23);

//  cout << "creo es esto: " << aval << "   " << bval << endl;

}

void Integration::stochasticTauIntegration(double &aval, double &bval, double tau) {

  int v1;
  double o13, o14, o24, o23, v13, v14, v23, v24, factor;
  o13 = 0.0;
  o14 = 0.0;
  o23 = 0.0;
  o24 = 0.0;
  v13 = 0.0;
  v14 = 0.0;
  v23 = 0.0;
  v24 = 0.0;

  for (int u=0; u<readData.numberOfOccOrb; u++) {
    factor = exp(tau*readData.Occ[u]);
 //   cout << "voy con:   " << readData.Occ[u] << endl;
    o13 += MosWP1[u] * MosWP3[u] * factor;
    o14 += MosWP1[u] * MosWP4[u] * factor;
    o23 += MosWP2[u] * MosWP3[u] * factor;
    o24 += MosWP2[u] * MosWP4[u] * factor;
  }

  v1 = readData.numberOfOccOrb;

//  cout << "Aqui hay cambio:   " << v1 << endl;
  for (int v=0; v<readData.numberOfVirOrb; v++) {
    factor = exp(-tau*readData.Vir[v]);
//    cout << "voy con:   " << readData.Vir[v] << endl;
    v13 += MosWP1[v1] * MosWP3[v1] * factor;
    v14 += MosWP1[v1] * MosWP4[v1] * factor;
    v23 += MosWP2[v1] * MosWP3[v1] * factor;
    v24 += MosWP2[v1] * MosWP4[v1] * factor;
    v1++;
  }
  aval = o13*o24*v13*v24 + o14*o23*v14*v23;
  bval = o14*o23*v13*v24 + o13*o24*v14*v23;


//  if (aval >= 1000 || bval >= 1000) 
//    cout << "CuidadOOOOOOO:  " << aval << "   " << bval << endl;

}


void Integration::generateMois(vector<double>& mois, double x, double y, double z) {
  double rr;
  double distx, disty, distz;
  double facx, facy, facz, expo, coef;
  int ang[3];
  int angMom;
  int funcType;
  int k, l;
  vector<double> gfs;
  gfs.resize(readData.numberOfMolOrb);

	l = 0;
	for (int i=0; i<readData.lastBasis.size(); i++) { 
		k = 0;
		for (int j=0; j<readData.lastBasis[i].shellSize.size(); j++) {
			gfs[l] = 0.0;
			for (int n=0; n<readData.lastBasis[i].shellSize[j]; n++) {
      	distx = x - readData.lastBasis[i].nucleusCoord[0];
      	disty = y - readData.lastBasis[i].nucleusCoord[1];
      	distz = z - readData.lastBasis[i].nucleusCoord[2];
      	rr = distx*distx + disty*disty + distz*distz;
				expo = readData.lastBasis[i].basisExponents[k]*rr;
				coef = readData.lastBasis[i].basisCoefficients[k];
				expo = exp(-expo);
				angMom = readData.lastBasis[i].angularMoment[k];
    		funcType = readData.lastBasis[i].functionType[k];
    		basis.angularMoment(ang, angMom, funcType);
    		facx = pow(distx, ang[0]);
    		facy = pow(disty, ang[1]);
    		facz = pow(distz, ang[2]);
				gfs[l] += facx*facy*facz*expo*coef;
				k++;
			}
			l++;
		}
	}

	for (int i=0; i<readData.numberOfMolOrb; i++) {
		k = 0;
		mois[i] = 0.0;
		for (int j=0; j<gfs.size(); j++) 
			mois[i] += gfs[j] * coeff.molecularOrbitalData[i].orbitalCoeff[j];
	}
/*
	for (int i=0; i<readData.lastBasis.size(); i++)
    for (int j=0; j<readData.lastBasis[i].basisExponents.size(); j++) {
      distx = x - readData.lastBasis[i].nucleusCoord[0];
      disty = y - readData.lastBasis[i].nucleusCoord[1];
      distz = z - readData.lastBasis[i].nucleusCoord[2];
      rr = distx*distx + disty*disty + distz*distz;
			expo = readData.lastBasis[i].basisExponents[j]*rr;
			if (expo > 45.0) 
				gfs[k] = 0.0;
			else {
    		expo = exp(-expo);
    		angMom = readData.lastBasis[i].angularMoment[j];
    		funcType = readData.lastBasis[i].functionType[j];
    		basis.angularMoment(ang, angMom, funcType);
    		facx = pow(distx, ang[0]);
    		facy = pow(disty, ang[1]);
    		facz = pow(distz, ang[2]);
				cout << "pin:  " << k << "   " << ang[0] << "  " << ang[1] << "  " << ang[2] << endl;
    		gfs[k] = facx*facy*facz*expo*readData.lastBasis[i].basisCoefficients[j];
			}
			k++;
    }

  for (int i=0; i<readData.numberOfMolOrb; i++) {
    k = 0;
		mois[i] = 0.0;
    for (int h=0; h<readData.lastBasis.size(); h++) 
      for (int m=0; m<readData.lastBasis[h].shellSize.size(); m++) 
        for (int n=0; n<readData.lastBasis[h].shellSize[m]; n++) {
					cout << "naa: " << k  << "   " << gfs[k] << "  " << readData.lastBasis[h].lastBasisCoeff[i].bC[m][n] << endl;
          mois[i] += gfs[k] * readData.lastBasis[h].lastBasisCoeff[i].bC[m][n];
          k++;
        }
				cout << "MIERDA: " << i << "   " << mois[i] << endl;
      }
*/
}

void Integration::generateMoisSpherical (vector<double>& gfs, double x, double y, double z) {
  double rr;
  double distx, disty, distz;
  double fac, expo, coef;
  int ang[3];
  int angMom;
  int funcType;
  int k, l;


  l = 0;
  for (int i=0; i<readData.lastBasis.size(); i++) {
    k = 0;
    distx = x - readData.lastBasis[i].nucleusCoord[0];
    disty = y - readData.lastBasis[i].nucleusCoord[1];
    distz = z - readData.lastBasis[i].nucleusCoord[2];
    rr = distx*distx + disty*disty + distz*distz;
//    cout << "ya pues : " << rr << "   " << distx << "   " << disty << "   " << distz << endl;
    for (int j=0; j<readData.lastBasis[i].shellSize.size(); j++) {
      gfs[l] = 0.0;
      for (int n=0; n<readData.lastBasis[i].shellSize[j]; n++) {
        expo = readData.lastBasis[i].basisExponents[k]*rr;
        coef = readData.lastBasis[i].basisCoefficients[k];
        expo = exp(-expo);
			  angMom = readData.lastBasis[i].angularMoment[k];
    		funcType = readData.lastBasis[i].functionType[k];
        fac = SphericalEvaluation(distx, disty, distz, rr, angMom, funcType);
        gfs[l] += coef*fac*expo;
//        cout << "ah ya: " << l << "  " << k << "   "  << gfs[l] << endl;
        k++;
      }
   //   cout << "chucha: " << l << "   " << gfs[l] << endl;
      l++;
    }
  }
  
/*
  for (int i=0; i<readData.numberOfMolOrb; i++) {
    k = 0;
    mois[i] = 0.0;
    for (int j=0; j<gfs.size(); j++) {
      mois[i] += gfs[j] * coeff.molecularOrbitalData[i].orbitalCoeff[j];
    }
  }
*/
}


double Integration::SphericalEvaluation(double x, double y, double z, double r, int ang, int type) {

  double result;
  if (ang == 0) {
    result = 1.0;
  } else if  (ang == 1) {
    switch (type) {
      case 0 : result = x;
               break;
      case 1 : result = y;
               break;
      case 2 : result = z;
               break;
      default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 2) {
    switch (type) {
      case 0 : result = sc[0]*x*y;
               break;
      case 1 : result = sc[0]*y*z;
               break;
      case 2 : result = 0.5*(2.0*z*z - x*x - y*y);
               break;
      case 3 : result = -sc[0]*x*z;
               break;
      case 4 : result = 0.5*sc[0]*(x*x-y*y);
               break;
      default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 3) {
    switch (type) {
      case 0: result = 0.5 * sc[1] * (3.0*x*x - y*y)*y;
              break;
      case 1: result = sc[2] * x * y * z;
              break;
      case 2: result = 0.5 * sc[3] * (4.0*z*z - x*x - y*y)*y;
              break;
      case 3: result = 0.5 * (2.0*z*z - 3.0 *(x*x + y*y))*z;
              break;
      case 4: result = 0.5 * sc[3] * (4.0*z*z - x*x - y*y)*x;
              break;
      case 5: result = 0.5 * sc[2] * (x*x - y*y)*z;
              break;
      case 6: result = 0.5 * sc[1] * (x*x - 3.0*y*y)*x;
              break;
      default: cout << "Check the type of primitive." << endl;
    }  
  } else if (ang == 4) {
    switch (type) {
      case 0: result = 0.5 * sc[4] * (x*x*x*y - y*y*x*y);
               break;
      case 1: result = 0.5 * sc[5] * (3.0*x*x*y*z - y*y*y*z);
              break;
      case 2: result = 0.5 * sc[6] * (6.0*z*z - y*y - x*x)*x*y;
              break;
      case 3: result = 0.5 * sc[1] * (4.0*z*z - 3.0*(x*x + y*y))*y*z;
              break;
      case 4: result = 4.375*z*z*z*z - 3.75*z*z*r*r + 0.375*r*r*r*r;
              break;
      case 5: result = 0.5 * sc[1] * (4.0*z*z - 3.0*(x*x + y*y))*x*z;
              break;
      case 6: result = 0.25 * sc[6] * (6.0*z*z - y*y - x*x)*(x*x - y*y);
              break;
      case 7: result = 0.5 * sc[5] * (x*x - 3.0*y*y)*x*z;
              break;
      case 8: result = 0.125 * sc[4] * (x*x*x*x - 6.0*x*x*y*y + y*y*y*y);
              break;
    default: cout << "Check the type of primitive." << endl;
    }
  }

  return result; 
}


void Integration::gaussKronrodEval() {

	int a1;
	double tau1, tau2, fac;

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

	for (int t=0; t<21; t++) {
		tauExpos[t].occExpo.resize(readData.numberOfOccOrb);
		tauExpos[t].virExpo.resize(readData.numberOfVirOrb);
		for (int i=0; i<readData.numberOfOccOrb; i++) { 
			tauExpos[t].occExpo[i] = exp(gk.tp[t]*readData.Occ[i]);
		}
		a1 = readData.numberOfOccOrb;
		for (int a=0; a<readData.numberOfVirOrb; a++) {
			tauExpos[t].virExpo[a] = exp(-gk.tp[t]*readData.Vir[a]);
			a1++;
		}
		gk.denom[t] = 1.0/(gk.tau[t]*gk.tau[t]);		
	}

/*	
	double result, result1;
	result1 = 0.0;
	result = 1.0;
	for (int j=0; j<21; j++) {
		result = 1.0;
		for (int i=0; i<readData.numberOfOccOrb; i++) 
			result *= tauExpos[j].occExpo[i];
			result1 += gk.wj[j]/(gk.tau[j]*gk.tau[j]) * result;
	}
	result1 *= 0.5;

	for (int j=0; j<readData.numberOfOccOrb; j++) {
  	fac += coeff.molecularOrbitalData[j].orbitalEnergie;
 		cout << coeff.molecularOrbitalData[j].orbitalEnergie << endl;
		}

	for (int j=0; j<20; j++) {
		tauExpos[j].occExpo.resize(readData.numberOfOccOrb);
		tauExpos[j].virExpo.resize(readData.numberOfVirOrb);
		for (int i=0; i<readData.numberOfOccOrb; i++) {
			tau1 = coeff.molecularOrbitalData[i].orbitalEnergie;
			tauExpos[j].occExpo[i] = -lg.wj[j] / tau1;
		}
		a1= readData.numberOfOccOrb;
		for (int a=0; a<readData.numberOfVirOrb; a++) {
			tau2 = coeff.molecularOrbitalData[a1].orbitalEnergie;
			tauExpos[j].virExpo[a] = lg.wj[j] / tau2;
			a1++;
		}
	}
	
	double integral = 0.0;

	for (int i=0; i<20; i++)
		integral += lg.wj[i];
	integral /= (-fac);

	cout << "the result is: " << setprecision(12) << integral << "    " << result1 << endl;
*/
}


