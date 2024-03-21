#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <algorithm>
#include "Spinors.h"
using namespace std;

double relBasis::complexReal(double a, double b, double c, double d) {

	double ac, bd, result;
	
	ac = a * c;
	bd = b * d;
  result = ac - bd;

	return result;
}

double relBasis::complexIm(double a, double b, double c, double d) {

  double ad, bc, result;

  ad = a * d;
  bc = b * c;
  result = ad + bc;

  return result;
}

void relBasis::source() {

	for (int j=0; j<4; j++)
  	for (int i=0; i<4; i++) {
    	bI.moSpi[j].Occvalmo[i].valuesMos[0].resize(bI.occOrb);
    	bI.moSpi[j].Occvalmo[i].valuesMos[1].resize(bI.occOrb);
    	bI.moSpi[j].Virvalmo[i].valuesMos[0].resize(bI.virOrb);
    	bI.moSpi[j].Virvalmo[i].valuesMos[1].resize(bI.virOrb);
		}


}

void relBasis::mofunc(double xinp, double yinp, double zinp, int spi) {

	double x, y, z, ar;
	int ibx, il, hi;
	int p1, p2, p3, pmn, pmm;
	double rsqrd, dcmpxxya1, dcmpxxyb1, dcmpxxya0, dcmpxxyb0;
	double texp, texp1, texp2, cc[4], xsign;
	
	for (int i=0; i<4; i++) 
		for (int j=0; j<bI.rdim; j++) {
			bI.valao[i].val[0][j] = 0.0;
			bI.valao[i].val[1][j] = 0.0;
		}


	for (int iatm=0; iatm<rG.atmNum; iatm++) {
		ibx = rG.geoRel[iatm].atomicNum;

		x = xinp - rG.geoRel[iatm].nucCoor[0];
		y = yinp - rG.geoRel[iatm].nucCoor[1];
		z = zinp - rG.geoRel[iatm].nucCoor[2];
			
		rsqrd = x*x + y*y + z*z;
	
		dcmpxxya0 = x;
		dcmpxxya1 = y;
		
		dcmpxxyb0 = x;
		dcmpxxyb1 = -y;

		bI.sph[0][0] = 1.0;
		bI.sph[1][0] = 0.0;
		
		bI.sph[0][3] = dcmpxxya0 / bI.sqV[1];
		bI.sph[1][3] = dcmpxxya1 / bI.sqV[1];
	
		bI.sph[0][2] = z;
		bI.sph[1][2] = 0.0;
		
		bI.sph[0][1] = -bI.sph[0][3];
		bI.sph[1][1] = bI.sph[1][3];	
			
		for (int iangl=0; iangl<=bI.bD[ibx-1].l; iangl++) {
	
			p1 =  iangl   * iangl    + 1;
     	p2 = (iangl-1)*(iangl-1) + 1;
     	p3 = (iangl+1)*(iangl+1) + 1;
	
			if(iangl>0) 
				make_sph(iangl, dcmpxxya0, dcmpxxya1, dcmpxxyb0, dcmpxxyb1, z);

			for (int ik=0; ik<bI.bD[ibx-1].cg[iangl].k; ik++) {

				ar = bI.bD[ibx-1].cg[iangl].ak[ik] * rsqrd;
				texp = exp(-ar);
				texp1 = -2.0*texp*bI.bD[ibx-1].cg[iangl].ak[ik];
				texp2 = texp*(2*iangl+1) + texp1 * rsqrd;
				il = bI.bD[ibx-1].cg[iangl].n;
				for (int in=0; in<il; in++) {
					if (iangl != 0) {
						cc[0] = bI.bD[ibx-1].cg[iangl].RL1[ik*il+in] * texp;
						cc[1] = bI.bD[ibx-1].cg[iangl].RS1[ik*il+in] * texp2;
						cc[2] = bI.bD[ibx-1].cg[iangl].RL2[ik*il+in] * texp;
						cc[3] = bI.bD[ibx-1].cg[iangl].RS2[ik*il+in] * texp1;
					} else {
						cc[0] = 0.0;
						cc[1] = 0.0;
						cc[2] = bI.bD[ibx-1].cg[iangl].RL1[ik*il+in] * texp;
						cc[3] = bI.bD[ibx-1].cg[iangl].RS1[ik*il+in] * texp1;
					}	
					
					for (int im=0; im<=2*iangl+1; im++) {

						if (2*iangl>im) {
							pmn = rG.geoRel[iatm].rpsb[im+2*iangl*iangl] + in;
							if (iangl-im-1 > 0)	
								xsign = pow(-1.0,iangl-im-1);
							else
								xsign = 1.0;

							if (2*iangl >= im && im >= 0) {
								bI.valao[0].val[0][pmn-1] += bI.cG[3].cgc[im][iangl]*bI.sph[0][p1+im-1]*cc[0]*xsign/bI.sqV[2*iangl];
								bI.valao[0].val[1][pmn-1] += bI.cG[3].cgc[im][iangl]*bI.sph[1][p1+im-1]*cc[0]*xsign/bI.sqV[2*iangl];
							}

							if (2*iangl >= im+1 && im+1>=0) {
								bI.valao[1].val[0][pmn-1] += bI.cG[2].cgc[im][iangl]*bI.sph[0][p1+im]*cc[0]*xsign/bI.sqV[2*iangl]; 
								bI.valao[1].val[1][pmn-1] += bI.cG[2].cgc[im][iangl]*bI.sph[1][p1+im]*cc[0]*xsign/bI.sqV[2*iangl]; 
							}
							if (2*(iangl-1) >= im-1 && im-1 >= 0) {
								bI.valao[2].val[0][pmn-1] += bI.cG[0].cgc[im][iangl-1]*bI.sph[0][p2+im-2]*cc[1]*xsign/bI.sqV[2*iangl];
								bI.valao[2].val[1][pmn-1] += bI.cG[0].cgc[im][iangl-1]*bI.sph[1][p2+im-2]*cc[1]*xsign/bI.sqV[2*iangl];
							}
							if (2*(iangl-1) >= im && im >= 0) {
								bI.valao[3].val[0][pmn-1] -= bI.cG[1].cgc[im][iangl-1]*bI.sph[0][p2+im-1]*cc[1]*xsign/bI.sqV[2*iangl]; 
								bI.valao[3].val[1][pmn-1] -= bI.cG[1].cgc[im][iangl-1]*bI.sph[1][p2+im-1]*cc[1]*xsign/bI.sqV[2*iangl]; 
							}		
						} // first if 2*iangl > im

						pmm = rG.geoRel[iatm].rpsb[im+ 2*iangl+ 2*iangl*iangl] + in;
						
						if (iangl-im+1 > 0)	
							xsign = pow(-1.0,iangl-im+1);
						else
							xsign = 1.0;

	// LARGE-----------------------
	// La
 	//C+j M * M-1/2
						if (2*iangl >= im-1 && im-1 >= 0) {
							bI.valao[0].val[0][pmm-1] += bI.cG[0].cgc[im][iangl]*bI.sph[0][p1+im-2]*cc[2]*xsign/bI.sqV[2*iangl];
							bI.valao[0].val[1][pmm-1] += bI.cG[0].cgc[im][iangl]*bI.sph[1][p1+im-2]*cc[2]*xsign/bI.sqV[2*iangl];
						}
						if (2*iangl >= im && im >= 0) {
							bI.valao[1].val[0][pmm-1] -= bI.cG[1].cgc[im][iangl]*bI.sph[0][p1+im-1]*cc[2]*xsign/bI.sqV[2*iangl];
							bI.valao[1].val[1][pmm-1] -= bI.cG[1].cgc[im][iangl]*bI.sph[1][p1+im-1]*cc[2]*xsign/bI.sqV[2*iangl];
						}
						if (2*(iangl+1) >= im && im >= 0) {
							bI.valao[2].val[0][pmm-1] += bI.cG[3].cgc[im][iangl+1]*bI.sph[0][p3+im-1]*cc[3]*xsign/bI.sqV[2*iangl];
							bI.valao[2].val[1][pmm-1] += bI.cG[3].cgc[im][iangl+1]*bI.sph[1][p3+im-1]*cc[3]*xsign/bI.sqV[2*iangl];
			//					cout << "pmm = " << pmm  << "  im =  " << im << " iangl+1  =  " << iangl+1 << "            " <<  bI.cG[3].cgc[im][iangl+1] << "            " << bI.valao[2].val[0][pmm-1] << "      " << bI.valao[2].val[1][pmm-1] << endl;
						}
						if (2*(iangl+1) >= im+1 && im+1 >= 0) {
							bI.valao[3].val[0][pmm-1]  += bI.cG[2].cgc[im][iangl+1]*(bI.sph[0][p3+im]*cc[3]*xsign/bI.sqV[2*iangl]);
							bI.valao[3].val[1][pmm-1]  += bI.cG[2].cgc[im][iangl+1]*bI.sph[1][p3+im]*cc[3]*xsign/bI.sqV[2*iangl];
						}
					} // for im	
				} // for in	
			} // for ik
		} // for iangl
	} // for iatm	
	
	////--------------------------sigue lo del valmo

  for (int j=0; j<4; j++) {
  	for (int i=0; i<bI.occOrb; i++) {
      bI.moSpi[spi].Occvalmo[j].valuesMos[0][i] = 0.0;
      bI.moSpi[spi].Occvalmo[j].valuesMos[1][i] = 0.0;
		}
  	for (int i=0; i<bI.virOrb; i++) {
      bI.moSpi[spi].Virvalmo[j].valuesMos[0][i] = 0.0;
      bI.moSpi[spi].Virvalmo[j].valuesMos[1][i] = 0.0;
		}
	}

	generateSpinors(spi);

}

void relBasis::generateSpinors (int si) {


	int hi;
	double totalR, totalI;

	hi = 0;
	for (int im=bI.rdim; im<bI.rdim+bI.occOrb; im++) { 
//	for (int im=bI.rdim; im<bI.rdim+10; im++) { 
		for (int in=0; in<bI.rdim; in++) { 

      // This part makes the real part
		  bI.moSpi[si].Occvalmo[0].valuesMos[0][hi] += complexReal(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
      bI.moSpi[si].Occvalmo[1].valuesMos[0][hi] += complexReal(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);;
      bI.moSpi[si].Occvalmo[2].valuesMos[0][hi] += complexReal(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
			bI.moSpi[si].Occvalmo[3].valuesMos[0][hi] += complexReal(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
      // This part makes the imaginary part
      bI.moSpi[si].Occvalmo[0].valuesMos[1][hi] += complexIm(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
			bI.moSpi[si].Occvalmo[1].valuesMos[1][hi] += complexIm(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
			bI.moSpi[si].Occvalmo[2].valuesMos[1][hi] += complexIm(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
 			bI.moSpi[si].Occvalmo[3].valuesMos[1][hi] += complexIm(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
    }
	/*	cout << bI.moSpi[si].Occvalmo[0].valuesMos[0][hi] << "   " << bI.moSpi[si].Occvalmo[0].valuesMos[1][hi] << endl; 
		cout << bI.moSpi[si].Occvalmo[1].valuesMos[0][hi] << "   " << bI.moSpi[si].Occvalmo[1].valuesMos[1][hi] << endl; 
		cout << "TOTAL: " << bI.moSpi[si].Occvalmo[2].valuesMos[0][hi] << "   " << bI.moSpi[si].Occvalmo[2].valuesMos[1][hi] << endl; 
		cout << bI.moSpi[si].Occvalmo[3].valuesMos[0][hi] << "   " << bI.moSpi[si].Occvalmo[3].valuesMos[1][hi] << endl; 
*/
//		totalR = bI.moSpi[si].Occvalmo[0].valuesMos[0][hi] + bI.moSpi[si].Occvalmo[1].valuesMos[0][hi] + bI.moSpi[si].Occvalmo[2].valuesMos[0][hi] + bI.moSpi[si].Occvalmo[3].valuesMos[0][hi];
//		totalI = bI.moSpi[si].Occvalmo[0].valuesMos[1][hi] + bI.moSpi[si].Occvalmo[1].valuesMos[1][hi] + bI.moSpi[si].Occvalmo[2].valuesMos[1][hi] + bI.moSpi[si].Occvalmo[3].valuesMos[1][hi]; 
//		cout << im << "   " << totalR << "   " << totalI << endl;
	hi++;
	} 

  hi = 0;
  for (int im=bI.rdim+bI.occOrb; im<bI.rdim*2; im++) {
    for (int in=0; in<bI.rdim; in++) {

      // This part makes the real part
      bI.moSpi[si].Virvalmo[0].valuesMos[0][hi] += complexReal(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
      bI.moSpi[si].Virvalmo[1].valuesMos[0][hi] += complexReal(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);;
      bI.moSpi[si].Virvalmo[2].valuesMos[0][hi] += complexReal(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
      bI.moSpi[si].Virvalmo[3].valuesMos[0][hi] += complexReal(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);

      // This part makes the imaginary part
      bI.moSpi[si].Virvalmo[0].valuesMos[1][hi] += complexIm(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
      bI.moSpi[si].Virvalmo[1].valuesMos[1][hi] += complexIm(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
      bI.moSpi[si].Virvalmo[2].valuesMos[1][hi] += complexIm(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
      bI.moSpi[si].Virvalmo[3].valuesMos[1][hi] += complexIm(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
    }
  hi++;
  }

}


void relBasis::make_sph(int mxl, double xa0, double xa1, double xb0, double xb1, double z) {

	int il, im, p0, p1, p2;

	il = mxl + 1;
	p0 = (il+1)*(il+1)+1;
  p1 = il*(il-1)+1;
	
	p0--;

	bI.sph[0][p0-1] = bI.sqV[2*il-2]/bI.sqV[2*il-1] * complexReal(xa0, xa1, bI.sph[0][p1+il-2], bI.sph[1][p1+il-2]);
	bI.sph[1][p0-1] = bI.sqV[2*il-2]/bI.sqV[2*il-1] * complexIm(xa0, xa1, bI.sph[0][p1+il-2], bI.sph[1][p1+il-2]);
	
  im = il-1;

	for (int i=im; i>=1; i--) {
		p0--;
		if ((il-i-1) != 0) { 
			bI.sph[0][p0-1] = (bI.sqV[il+i-1]/bI.sqV[il-i-1]) * z * bI.sph[0][p1+i-1] - (bI.sqV[il-i-2]/bI.sqV[il-i-1]) * complexReal(xb0, xb1, bI.sph[0][p1+i], bI.sph[1][p1+i]);
			bI.sph[1][p0-1] = (bI.sqV[il+i-1]/bI.sqV[il-i-1]) * z * bI.sph[1][p1+i-1] - (bI.sqV[il-i-2]/bI.sqV[il-i-1]) * complexIm(xb0, xb1, bI.sph[0][p1+i], bI.sph[1][p1+i]);
		} else {
			bI.sph[0][p0-1] = bI.sqV[il+i-1]/bI.sqV[il-i-1] * z * bI.sph[0][p1+i-1];
			bI.sph[1][p0-1] = bI.sqV[il+i-1]/bI.sqV[il-i-1] * z * bI.sph[1][p1+i-1];
		}
	}

	p0--;
	
	bI.sph[0][p0-1] = z * bI.sph[0][p1-1] - (bI.sqV[il-2]/bI.sqV[il-1]) * complexReal(xb0, xb1, bI.sph[0][p1], bI.sph[1][p1]);
	bI.sph[1][p0-1] = z * bI.sph[1][p1-1] - (bI.sqV[il-2]/bI.sqV[il-1]) * complexIm(xb0, xb1, bI.sph[0][p1], bI.sph[1][p1]);

	p2 = p0;


	for (int i=0; i<il; i++) {
		p0--;
		p2++;
		if (i%2 == 0) {
			bI.sph[0][p0-1] = -bI.sph[0][p2-1];
			bI.sph[1][p0-1] = bI.sph[1][p2-1];
		} else { 
			bI.sph[0][p0-1] = bI.sph[0][p2-1];
			bI.sph[1][p0-1] = -bI.sph[1][p2-1];
		}
	}

}

void relBasis::readBasisRel(fstream& fileInfo, fstream& fileBasis, fstream& xyzFile, fstream& coeffsFile, fstream& orbEneFile) {
  
  rG.geomLecture(xyzFile);
  string file1 = "../h2.MP2OPT/H2.MP2.R4DInfo";
  string file2 = "../h2.MP2OPT/H2.MP2.BasisUT";
//  string file1 = "../H2O.MP2/h2o.MP2.R4DInfo";
//  string file2 = "../H2O.MP2/h2o.MP2.BasisUT";
//  string file1 = "../cuh.MP2OPT/cuh.MP2OPT.R4DInfo";
//  string file2 = "../cuh.MP2OPT/cuh.MP2OPT.BasisUT";
//  string file1 = "../cu2.MP2OPT/cu2.MP2.R4DInfo";
//  string file2 = "../cu2.MP2OPT/cu2.MP2.BasisUT";
//  string file1 = "../agh.MP2OPT/agh.MP2OTP.R4DInfo";
//  string file2 = "../agh.MP2OPT/agh.MP2OTP.BasisUT";
//  string file1 = "../ag2.MP2OPT/ag2.MP2OPT.R4DInfo";
//  string file2 = "../ag2.MP2OPT/ag2.MP2OPT.BasisUT";
//  string file1 = "../auh.MP2/auh.MP2.R4DInfo";
//  string file2 = "../auh.MP2/auh.MP2.BasisUT";
//  string file1 = "../au2.MP2/au2.MP2.R4DInfo";
//  string file2 = "../au2.MP2/au2.MP2.BasisUT";


  string RELA;
  int l, k, n, im, u, ik, j;

  ifstream infoFile(file1.c_str(), ios::in);

  if (!infoFile) {
    cerr << "The file R4DInfo can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  ifstream basisFile(file2.c_str(), ios::in);

  if (!basisFile) {
    cerr << "The file BasisUT can't be open." << endl;
    exit(EXIT_FAILURE);
  }
  
  infoFile >> bI.maxBid >> bI.maxl;
	
  infoFile.close();
  
//  cout << "max bid = " << bI.maxBid << endl;
//  cout << "max l = "  << bI.maxl << endl;

  bI.bD.resize(bI.maxBid);
  
  for (int j=0; j<bI.maxBid; j++) {
//    cout << "bid = " << j << endl;
    getline(basisFile, RELA);
//    cout << RELA << endl;
    if ( j != 0) {     
      getline(basisFile, RELA);
 //     cout << RELA << endl;
    }
    basisFile >> l;
    bI.bD[j].l = l - 1;
    bI.bD[j].cg.resize(bI.bD[j].l+1);
		for (int i=0; i<=bI.bD[j].l; i++) {
      basisFile >> k;
      basisFile >> n;
      basisFile >> im;
      bI.bD[j].cg[i].k = k;
      bI.bD[j].cg[i].n = n;
      bI.bD[j].cg[i].ak.resize(k);
      bI.bD[j].cg[i].RL1.resize(k*n);
      bI.bD[j].cg[i].RS1.resize(k*n);
      bI.bD[j].cg[i].RL2.resize(k*n);
      bI.bD[j].cg[i].RS2.resize(k*n);
      for (int ik=0; ik<k; ik++) {
        basisFile >> bI.bD[j].cg[i].ak[ik];
        for (int in=0; in<n; in++) { 
          basisFile >> bI.bD[j].cg[i].RL1[ik*n + in];
        }
      }
      for (int ik=0; ik<k; ik++) 
        for (int in=0; in<n; in++) {
          basisFile >> bI.bD[j].cg[i].RS1[ik*n + in];
        } 
      if (i != 0) {
        for (int ik=0; ik<k; ik++) 
          for (int in=0; in<n; in++) {
            basisFile >> bI.bD[j].cg[i].RL2[ik*n + in];
          }
        for (int ik=0; ik<k; ik++) 
          for (int in=0; in<n; in++) {
            basisFile >> bI.bD[j].cg[i].RS2[ik*n + in];
          }
        }
    }
  }
  
  basisFile.close();
  
  size_t num;
  num = 2*(10+bI.maxl);
  bI.sqV.resize(num);
  bI.sph[0].resize((bI.maxl+1)*(bI.maxl+1)+10);
  bI.sph[1].resize((bI.maxl+1)*(bI.maxl+1)+10);

  for (int i=0; i<num; i++) 
    bI.sqV[i] = sqrt(i+1);

  int mxBl, mxBl2;
  mxBl = bI.maxl + 1;
  mxBl2 = 2*mxBl + 1;

	for (int i=0; i<4; i++) {
		bI.cG[i].cgc.resize(mxBl2+1);
		for (int j=0; j<=mxBl2; j++)
			bI.cG[i].cgc[j].resize(mxBl+1);	
	}


  for (int i=0; i<=mxBl; i++) 
    for (int m=0; m<=2*i+1; m++) {
      bI.cG[0].cgc[m][i] = sqrt(max(0.0,m*1.0));
      bI.cG[1].cgc[m][i] = sqrt(max(0.0,1.0*(2*i-m+1)));
      bI.cG[2].cgc[m][i] = sqrt(max(0.0,1.0*(m+1)));
      bI.cG[3].cgc[m][i] = sqrt(max(0.0,1.0*(2*i-m)));
	}
   

  int mi;
	size_t num1;

  ik = 1;
  for (int i=0; i<rG.atmNum; i++) {
    j = rG.geoRel[i].atomicNum - 1;
    l = bI.bD[j].l;
    num1 = 2*((l+1)*(l+1))-1;
    rG.geoRel[i].rpsb.resize(num1, 0);
    im = 0;
    for (int m=0; m<2*l+1; m++) { 
      mi = int(m/2)*2 + 2;
      for (int n=0; n<mi; n++) {
        rG.geoRel[i].rpsb[im] = ik;
        ik += bI.bD[j].cg[(m+1)/2].n;
        im++;
      }    
    }
  }

  bI.rdim = ik - 1;

	for (int u=0; u<4; u++) 
		for (int l=0; l<2; l++) 
  		bI.valao[u].val[l].resize(bI.rdim);

  readOrbitalCoeff(coeffsFile);

  readOrbitalEne(orbEneFile);

  source();

//  cout << "rdim = " << bI.rdim << endl;

}

void relBasis::readOrbitalEne(fstream& fileOrbEne) {

	double enVal;
	int ignoreInt;

//	string file = "../H2O.MP2/h2o.MP2.r4dorbene";
	string file = "../h2.MP2OPT/H2.MP2.r4dorbene";
//    string file = "../cuh.MP2OPT/cuh.MP2OPT.r4dorbene";
//    string file = "../cu2.MP2OPT/cu2.MP2.r4dorbene";
//    string file = "../agh.MP2OPT/agh.MP2OTP.r4dorbene";
//    string file = "../ag2.MP2OPT/ag2.MP2OPT.r4dorbene";
//    string file = "../au2.MP2/au2.MP2.r4dorbene";
//    string file = "../auh.MP2/auh.MP2.r4dorbene";

//  ifstream eneFile(file.c_str(), ios::in);
  ifstream eneFile(file.c_str(), ios::binary);

  if (!eneFile) {
    cerr << "The file r4dorbene can't be open." << endl;
    exit(EXIT_FAILURE);
  }

	bI.orbEne.resize(bI.rdim);	

	eneFile.read((char*)&ignoreInt, 4);
	for (int i=0; i<bI.rdim; i++) 
		eneFile.read((char*)&enVal, 8);
	for (int i=0; i<bI.rdim; i++) {
		eneFile.read((char*)&enVal, 8);
		bI.orbEne[i] = enVal;
	}

	eneFile.close();

	ignoreInt = 0;
	enVal = 0.0;
/*
  	while (enVal <= 0.0) {
		enVal = bI.orbEne[ignoreInt];
		if (enVal < 0.0)
			ignoreInt++;
	};  
*/

	bI.occOrb = rG.nE;
	bI.virOrb = bI.rdim - bI.occOrb;

} 

void relGeometry::geomLecture(fstream& fileXYZ) {

	int j, t1, t2;
  string file = "../h2.MP2OPT/H2.MP2.AtomUT";
//  string file = "../H2O.MP2/h2o.MP2.AtomUT";
//  string file = "../cuh.MP2OPT/cuh.MP2OPT.AtomUT";
//  string file = "../cu2.MP2OPT/cu2.MP2.AtomUT";
//  string file = "../agh.MP2OPT/agh.MP2OTP.AtomUT";
//  string file = "../ag2.MP2OPT/ag2.MP2OPT.AtomUT";
//  string file = "../au2.MP2/au2.MP2.AtomUT";
//  string file = "../auh.MP2/auh.MP2.AtomUT";

  ifstream geomFile(file.c_str(), ios::in);

  if (!geomFile) {
    cerr << "The file AtomUT can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  geomFile >> atmNum;
  geoRel.resize(atmNum);
//  cout << "There are: " << atmNum << endl;
  nE = 0;
  for (int i=0; i<atmNum; i++) {
    geomFile >> geoRel[i].atomType >> geoRel[i].atomicNum >> geoRel[i].nucCoor[0] >> geoRel[i].nucCoor[1] >> geoRel[i].nucCoor[2] >> geoRel[i].charge;
    nE += int(geoRel[i].charge);
  } 
  geomFile.close();

	difAtoms.resize(1);
	difAtoms[0] = geoRel[0].charge;
	if (atmNum > 1) {
		for (int i=0; i<atmNum; i++) {
      t1 = geoRel[i].charge;
      j = 0;
      while (j < difAtoms.size()) {
        t2 = difAtoms[j];
        if (t1 == t2)
          j = atmNum;
        else
          j++;
      }
      if (j == difAtoms.size())
        difAtoms.push_back(t1);
    }
  }

  nDifAtoms = difAtoms.size();

}

void relBasis::readOrbitalCoeff(fstream& readFromFile) {

  int ignoreInt;
	double orb0, orb1, orb2, orb3;
	size_t fac;
	
	fac = 2*bI.rdim;
	bI.mo[0].moS.resize(fac);
  bI.mo[1].moS.resize(fac);

  for (int i=0; i<fac; i++) {
    bI.mo[0].moS[i].moRe.resize(fac);
    bI.mo[1].moS[i].moRe.resize(fac);
    bI.mo[0].moS[i].moIm.resize(fac);
    bI.mo[1].moS[i].moIm.resize(fac);
		for (int j=0; j<fac; j++) {
			bI.mo[0].moS[i].moRe[j] = 0.0;
			bI.mo[1].moS[i].moRe[j] = 0.0;
			bI.mo[0].moS[i].moIm[j] = 0.0;
			bI.mo[1].moS[i].moIm[j] = 0.0;
		}

  }

  string file = "../h2.MP2OPT/salidaH2MP2";
//  string file = "../H2O.MP2/salidaH2OMP2";
//  string file = "../cuh.MP2OPT/salidaCuHMP2";
//  string file = "../cu2.MP2OPT/salidaCu2MP2";
//    string file = "../agh.MP2OPT/salidaAgHMP2";
//    string file = "../ag2.MP2OPT/salidaAg2MP2";
//    string file = "../auh.MP2/salidaAuHMP2";
//    string file = "../au2.MP2/salidaAu2MP2";


  ifstream input;
  input.open(file.c_str(), ios::in);

  if (!input.is_open()) {
    cerr << "The file coeff can't be open." << endl;
    exit(EXIT_FAILURE);
  }
  
  string line;
  char character1, character2, character3;
  double value1, value0;
  for (int k=0; k<2; k++) 
    for (int i=0; i<fac; i++) 
      for (int j=0; j<fac; j++) {
        getline(input, line);
        istringstream lineCheck(line); 
        lineCheck >> character1 >> value0 >> character2 >> value1 >> character3;
        bI.mo[k].moS[i].moRe[j] = value0;
        bI.mo[k].moS[i].moIm[j] = value1;
  }

 input.close();

}


void relBasis::weightFunctionData(fstream& MCFile) {

	int j, atom;
	size_t l = rG.atmNum;
	size_t t = rG.nDifAtoms;
	wg.wV.resize(t);
  wg.wD.resize(l);

//  string file = "MCDataAgH";
  string file = "MCDataH2";
//  string file = "MCDataCuH";
//  string file = "MCDataAu2";

  ifstream MC(file.c_str(), ios::in);

  if (!MC) {
    cerr << "The MC file can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  for (int i=0; i<t; i++) {
    MC >> wg.wV[i].coeffWeight[0];
    MC >> wg.wV[i].coeffWeight[1];
    MC >> wg.wV[i].expWeight[0];
    MC >> wg.wV[i].expWeight[1];
    wg.wV[i].normWeight[0] = wg.normalize(wg.wV[i].expWeight[0]);
    wg.wV[i].normWeight[1] = wg.normalize(wg.wV[i].expWeight[1]);
 //   cout << "normalize: " << wg.wV[i].normWeight[0] *  wg.wV[i].coeffWeight[0]<< "  " << wg.wV[i].normWeight[1]*wg.wV[i].coeffWeight[1] << endl;
    wg.wV[i].atomType = rG.difAtoms[i];
  }

  MC >> nWalkers;
  MC >> nSteps;
  MC >> Nbsize;
  
  MC.close();

  for (int i=0; i<l; i++) {
    atom = rG.geoRel[i].charge;
    j = 0;
    while (j<rG.difAtoms.size()) {
      if (atom == wg.wV[j].atomType) {
        wg.wD[i].atomType = atom;
        copy(begin(rG.geoRel[i].nucCoor), end(rG.geoRel[i].nucCoor), begin(wg.wD[i].nucleusCoord));
        copy(begin(wg.wV[j].coeffWeight), end(wg.wV[j].coeffWeight), begin(wg.wD[i].coeffWeight));
        copy(begin(wg.wV[j].expWeight), end(wg.wV[j].expWeight), begin(wg.wD[i].expWeight));
        copy(begin(wg.wV[j].normWeight), end(wg.wV[j].normWeight), begin(wg.wD[i].normWeight));
        j = t;
      }
      else
        j++;
    }
  }


/*  for (int i=0; i<wg.wD.size(); i++) {
    cout << wg.wD[i].atomType << endl;
    cout << wg.wD[i].nucleusCoord[0] << " " << wg.wD[i].nucleusCoord[1] << " " << wg.wD[i].nucleusCoord[2] << endl;
    for (int j=0; j<wg.wD[i].coeffWeight.size(); j++)
      cout << " coeff: " << j << " " << wg.wD[i].coeffWeight[j] << " " << " expo: " << wg.wD[i].expWeight[j] << endl;
  }
*/
	
}


double relBasis::constant () {
  double result;
  double result1;
  double Ax, Bx, Ay, By, Az, Bz;
  double Rab;
  size_t t = wg.wD.size();

  result = 0.0;
  for (int i=0; i<t; i++)
    for (int j=0; j<t; j++) {
      Ax = wg.wD[i].nucleusCoord[0];
      Ay = wg.wD[i].nucleusCoord[1];
      Az = wg.wD[i].nucleusCoord[2];

      Bx = wg.wD[j].nucleusCoord[0];
      By = wg.wD[j].nucleusCoord[1];
      Bz = wg.wD[j].nucleusCoord[2];

      Rab = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
      Rab = sqrt(Rab);

      result1 = 0.0;
      for (int u=0; u<wg.wD[i].coeffWeight.size(); u++)
        for (int v=0; v<wg.wD[j].coeffWeight.size(); v++)
          result1 += wg.wD[i].normWeight[u]*wg.wD[j].normWeight[v]*wg.wD[i].coeffWeight[u]*wg.wD[j].coeffWeight[v]*wg.weightIntegral(wg.wD[i].expWeight[u], wg.wD[j].expWeight[v], Rab);

    result += result1;
  }


  return result;
}

double relBasis::weightFunction (double x1, double y1, double z1, double x2, double y2, double z2, double r12) {

  double Ax, Ay, Az, Bx, By, Bz;
  double result, coef1, expo1, coef2, expo2;
  double X1a, Y1a, Z1a, X2b, Y2b, Z2b, R1a, R2b, Rab;
  double gf1, gf2;
  size_t t = wg.wD.size();

  gf1 = 0.0;
  gf2 = 0.0;

//  for (int i=0; i<t; i++) 
    for (int j=0; j<t; j++) {
      Ax = wg.wD[j].nucleusCoord[0];
      Ay = wg.wD[j].nucleusCoord[1];
      Az = wg.wD[j].nucleusCoord[2];

      Bx = wg.wD[j].nucleusCoord[0];
      By = wg.wD[j].nucleusCoord[1];
      Bz = wg.wD[j].nucleusCoord[2];

      X1a = x1 - Ax;
      Y1a = y1 - Ay;
      Z1a = z1 - Az;
      R1a = X1a*X1a + Y1a*Y1a + Z1a*Z1a;

      X2b = x2 - Bx;
      Y2b = y2 - By;
      Z2b = z2 - Bz;
      R2b = X2b*X2b + Y2b*Y2b + Z2b*Z2b;

 //     for (int u=0; u<wD[i].coeffWeight.size(); u++) 
      for (int v=0; v<wg.wD[j].coeffWeight.size(); v++) {
        coef2 = wg.wD[j].coeffWeight[v] * wg.wD[j].normWeight[v];
        expo2 = wg.wD[j].expWeight[v];
        coef1 = wg.wD[j].coeffWeight[v] * wg.wD[j].normWeight[v];
        expo1 = wg.wD[j].expWeight[v];
        gf1 += coef1 * exp(-expo1*R1a);
        gf2 += coef2 * exp(-expo2*R2b);
        }
    }

  result = gf1 * gf2 / r12;


  return result;
}

