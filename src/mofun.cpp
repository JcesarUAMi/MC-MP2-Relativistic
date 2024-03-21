#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <algorithm>
using namespace std;

struct geomInfo {
  string atomType;
  int atomicNum;
  double charge;
  vector<int> rpsb;
  array<double, 3> nucCoor;
};

class relGeometry {
  public:
    int atmNum;
    vector<geomInfo> geoRel;
    void geomLecture (fstream&);
};

struct cgto {
  int n;
  int k;
  vector<double> ak;
  vector<double> RL1;
  vector<double> RS1;
  vector<double> RL2;
  vector<double> RS2;
};

struct basisData {
  int l;
  vector<cgto> cg;
};

struct cgcValues {
  vector<double> cgc0;
  vector<double> cgc1;
  vector<double> cgc2;
  vector<double> cgc3;
};

struct coeffValues {
	vector<double> moRe;
	vector<double> moIm;
};

struct spiMos {
  vector<coeffValues> moS;
	vector<vector<double>> ov;
};

struct valaoStruct {
  array<vector<double>, 2> val;
};

struct valmoValues {
	array<vector<double>, 2> valuesMos;
};

struct basisInfo {
  int maxBid;
  int maxl;
	int ibegin;
	int iend;
  int delta;
  int rdim;
	int occOrb, virOrb;
  vector<basisData> bD;
  cgcValues cgc;
  array<vector<double>, 2> sph;
  vector<double> sqV;
  array<spiMos, 2> mo;
	vector<double> orbEne;
  array<valaoStruct, 4> valao;
	array<valmoValues, 4> valmo;
};

class relBasis {
  public:
    basisInfo bI;
    relGeometry rG;
    void readBasisRel (fstream&, fstream&, fstream&, fstream&, fstream&);
    void readOrbitalCoeff(fstream&);
		void readOrbitalEne(fstream&);
		void mofunc(double, double, double);
		void make_sph(int, double, double, double, double, double);
		double complexReal(double, double, double, double);
		double complexIm(double, double, double, double);
		void source();
};

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

	double x, y, z, wg;
	int im, in; 
	double vIm0, vIm, vRe0, vRe;
	double n, l;
	
	x = 0.1;
	y = 0.2;
	z = 0.3;
	
	bI.ibegin = bI.rdim;
	bI.iend = bI.rdim + 2;
	
	bI.delta = bI.iend - bI.ibegin;

	bI.mo[0].ov.resize(bI.delta);
	bI.mo[1].ov.resize(bI.delta);

  for (int i=0; i<4; i++) {
    bI.valmo[i].valuesMos[0].resize(bI.delta);
    bI.valmo[i].valuesMos[1].resize(bI.delta);
	}
	
	mofunc(x, y, z);

	for (int j=0; j<bI.delta; j++) {
		bI.mo[0].ov[j].resize(bI.delta);
		bI.mo[1].ov[j].resize(bI.delta);
		for (int i=0; i<4; i++) {
			bI.valmo[i].valuesMos[0][j] = 0.0;
			bI.valmo[i].valuesMos[1][j] = 0.0;
		}
	}

	for (int u=0; u<bI.delta; u++) 
		for (int l=0; l<bI.delta; l++) {
			bI.mo[0].ov[u][l] = 0.0;
			bI.mo[1].ov[u][l] = 0.0;
	}
	
	n = 200;
 	l = 10;
	wg = 1.25E-4;

	for (int i=0; i<n; i++) 
		for (int j=0; j<n; j++)
			for (int k=0; k<n; k++) {
		    x = l*((i-n/2))/n + 0.05;
				y = l*((j-n/2))/n + 0.05;
				z = l*((k-n/2))/n + 0.05;

				cout << "CAMBIO: " << x << "  " << y << "  " << z << endl;
			
				mofunc(x, y, z);

				for (int iv=0; iv<bI.delta; iv++)
					for (int jv=0; jv<bI.delta; jv++) 
						for (int m=0; m<4; m++) {
							bI.mo[0].ov[jv][iv] += wg * complexReal(bI.valmo[m].valuesMos[0][jv], -bI.valmo[m].valuesMos[1][jv], bI.valmo[m].valuesMos[0][iv], bI.valmo[m].valuesMos[1][iv]);
							bI.mo[1].ov[jv][iv] += wg * complexIm(bI.valmo[m].valuesMos[0][jv], -bI.valmo[m].valuesMos[1][jv], bI.valmo[m].valuesMos[0][iv], bI.valmo[m].valuesMos[1][iv]);
						}
	    }

	for (int i=0; i<bI.delta; i++)
		for (int j=0; j<bI.delta; j++)
		cout << i << "  " << j << "   " << bI.mo[0].ov[i][j] << endl;

}

void relBasis::mofunc(double xinp, double yinp, double zinp) {

	double x, y, z, ar;
	int ibx, il;
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
								bI.valao[0].val[0][pmn-1] += bI.cgc.cgc3[(iangl+1)*iangl+im]*bI.sph[0][p1+im-1]*cc[0]*xsign/bI.sqV[2*iangl];
								bI.valao[0].val[1][pmn-1] += bI.cgc.cgc3[(iangl+1)*iangl+im]*bI.sph[1][p1+im-1]*cc[0]*xsign/bI.sqV[2*iangl];
							}

							if (2*iangl >= im+1 && im+1>=0) {
								bI.valao[1].val[0][pmn-1] += bI.cgc.cgc2[(iangl+1)*iangl+im]*bI.sph[0][p1+im]*cc[0]*xsign/bI.sqV[2*iangl]; 
								bI.valao[1].val[1][pmn-1] += bI.cgc.cgc2[(iangl+1)*iangl+im]*bI.sph[1][p1+im]*cc[0]*xsign/bI.sqV[2*iangl]; 
							}
							if (2*(iangl-1) >= im-1 && im-1 >= 0) {
								bI.valao[2].val[0][pmn-1] += bI.cgc.cgc0[iangl*(iangl-1)+im]*bI.sph[0][p2+im-2]*cc[1]*xsign/bI.sqV[2*iangl];
								bI.valao[2].val[1][pmn-1] += bI.cgc.cgc0[iangl*(iangl-1)+im]*bI.sph[1][p2+im-2]*cc[1]*xsign/bI.sqV[2*iangl];
							}
							if (2*(iangl-1) >= im && im >= 0) {
								bI.valao[3].val[0][pmn-1] -= bI.cgc.cgc1[iangl*(iangl-1)+im]*bI.sph[0][p2+im-1]*cc[1]*xsign/bI.sqV[2*iangl]; 
								bI.valao[3].val[1][pmn-1] -= bI.cgc.cgc1[iangl*(iangl-1)+im]*bI.sph[1][p2+im-1]*cc[1]*xsign/bI.sqV[2*iangl]; 
							}		
						} // first if 2*iangl > im

						pmm = rG.geoRel[iatm].rpsb[im+ 2*iangl+ 2*iangl*iangl] + in;
						
						if (iangl-im+1 > 0)	
							xsign = pow(-1.0,iangl-im+1);
						else
							xsign = 1.0;

					int index;	
	// LARGE-----------------------
	// La
 	//C+j M * M-1/2
						if (2*iangl >= im-1 && im-1 >= 0) {
							bI.valao[0].val[0][pmm-1] += bI.cgc.cgc0[(iangl+1)*iangl+im]*bI.sph[0][p1+im-2]*cc[2]*xsign/bI.sqV[2*iangl];
							bI.valao[0].val[1][pmm-1] += bI.cgc.cgc0[(iangl+1)*iangl+im]*bI.sph[1][p1+im-2]*cc[2]*xsign/bI.sqV[2*iangl];
						}
						if (2*iangl >= im && im >= 0) {
							bI.valao[1].val[0][pmm-1] -= bI.cgc.cgc1[(iangl+1)*iangl+im]*bI.sph[0][p1+im-1]*cc[2]*xsign/bI.sqV[2*iangl];
							bI.valao[1].val[1][pmm-1] -= bI.cgc.cgc1[(iangl+1)*iangl+im]*bI.sph[1][p1+im-1]*cc[2]*xsign/bI.sqV[2*iangl];
						}
						if (2*(iangl+1) >= im && im >= 0) {
           
              if (im == 0) {
                if (iangl == 0)
                  index = 2;
                else
                  index = 2*(iangl*iangl+2);
              } else {
                if (iangl == 0)
                  index = im+2;
                else
                  index = 2*(iangl*iangl+2)+im;
              }

							bI.valao[2].val[0][pmm-1] += bI.cgc.cgc3[index]*bI.sph[0][p3+im-1]*cc[3]*xsign/bI.sqV[2*iangl];
							bI.valao[2].val[1][pmm-1] += bI.cgc.cgc3[index]*bI.sph[1][p3+im-1]*cc[3]*xsign/bI.sqV[2*iangl];
						}
						if (2*(iangl+1) >= im+1 && im+1 >= 0) {
							bI.valao[3].val[0][pmm-1]  += bI.cgc.cgc2[(iangl+1)*iangl+im]*(bI.sph[0][p3+im]*cc[3]*xsign/bI.sqV[2*iangl]);
							bI.valao[3].val[1][pmm-1]  += bI.cgc.cgc2[(iangl+1)*iangl+im]*bI.sph[1][p3+im]*cc[3]*xsign/bI.sqV[2*iangl];
						}
					} // for im	
				} // for in	
			} // for ik
		} // for iangl
	} // for iatm	
	
	////--------------------------sigue lo del valmo


	int hi;

  for (int i=0; i<bI.delta; i++) 
    for (int j=0; j<4; j++) {
      bI.valmo[j].valuesMos[0][i] = 0.0;
      bI.valmo[j].valuesMos[1][i] = 0.0;
    }

  double result;
	hi = 0;
	for (int im=bI.ibegin; im<bI.iend; im++) { 
		for (int in=0; in<bI.rdim; in++) { 

      // This part makes the real part
		  bI.valmo[0].valuesMos[0][hi] += complexReal(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
      bI.valmo[1].valuesMos[0][hi] += complexReal(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);;
      bI.valmo[2].valuesMos[0][hi] += complexReal(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
			bI.valmo[3].valuesMos[0][hi] += complexReal(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);

      // This part makes the imaginary part
      bI.valmo[0].valuesMos[1][hi] += complexIm(bI.valao[0].val[0][in], bI.valao[0].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
			bI.valmo[1].valuesMos[1][hi] += complexIm(bI.valao[1].val[0][in], bI.valao[1].val[1][in], bI.mo[0].moS[in].moRe[im], bI.mo[0].moS[in].moIm[im]);
			bI.valmo[2].valuesMos[1][hi] += complexIm(bI.valao[2].val[0][in], bI.valao[2].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
 			bI.valmo[3].valuesMos[1][hi] += complexIm(bI.valao[3].val[0][in], bI.valao[3].val[1][in], bI.mo[0].moS[in+bI.rdim].moRe[im], bI.mo[0].moS[in+bI.rdim].moIm[im]);
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
//  string file1 = "../h2o.r4dscf.R4DInfo";
//  string file2 = "../h2o.r4dscf.BasisUT";
  string file1 = "../h2.r4dscf.R4DInfo";
  string file2 = "../h2.r4dscf.BasisUT";

  string RELA;
  int l, k, n, im, u, ik, j;

  ifstream infoFile(file1.c_str(), ios::in);

  if (!infoFile) {
    cerr << "The file can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  ifstream basisFile(file2.c_str(), ios::in);

  if (!basisFile) {
    cerr << "The file can't be open." << endl;
    exit(EXIT_FAILURE);
  }
  
  infoFile >> bI.maxBid >> bI.maxl;
	
  infoFile.close();
  
  cout << "max bid = " << bI.maxBid << endl;
  cout << "max l = "  << bI.maxl << endl;

  bI.bD.resize(bI.maxBid);
  
  for (int j=0; j<bI.maxBid; j++) {
    cout << "bid = " << j << endl;
    getline(basisFile, RELA);
    cout << RELA << endl;
    if ( j != 0) {     
      getline(basisFile, RELA);
      cout << RELA << endl;
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

	if (mxBl <= 1) {
  	bI.cgc.cgc0.resize((mxBl+1)*mxBl2);
  	bI.cgc.cgc1.resize((mxBl+1)*mxBl2);
  	bI.cgc.cgc2.resize((mxBl+1)*mxBl2);
  	bI.cgc.cgc3.resize((mxBl+1)*mxBl2);
	} else {
  	bI.cgc.cgc0.resize((mxBl2+1)*mxBl);
  	bI.cgc.cgc1.resize((mxBl2+1)*mxBl);
  	bI.cgc.cgc2.resize((mxBl2+1)*mxBl);
  	bI.cgc.cgc3.resize((mxBl2+1)*mxBl);
	}

  u = 0;
  for (int i=0; i<=mxBl; i++) 
    for (int m=0; m<=2*i+1; m++) {
      bI.cgc.cgc0[u] = sqrt(m);
      bI.cgc.cgc1[u] = sqrt(2*i-m+1);
      bI.cgc.cgc2[u] = sqrt(m+1);
      if (2*i < m) 
        bI.cgc.cgc3[u] = sqrt(0);
      else
        bI.cgc.cgc3[u] = sqrt(2*i-m);
      u++;
    }
   
	bI.cgc.cgc0.resize(u);
	bI.cgc.cgc1.resize(u);
	bI.cgc.cgc2.resize(u);
	bI.cgc.cgc3.resize(u);

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
  cout << "rdim = " << bI.rdim << endl;

	for (int u=0; u<4; u++) 
		for (int l=0; l<2; l++) 
  		bI.valao[u].val[l].resize(bI.rdim);

  readOrbitalCoeff(coeffsFile);

  readOrbitalEne(orbEneFile);

}

void relBasis::readOrbitalEne(fstream& fileOrbEne) {

	double enVal;
	int ignoreInt;

	string file = "../h2.r4dscf.r4dorbene";
//	string file = "../h2o.r4dscf.r4dorbene";
//  ifstream eneFile(file.c_str(), ios::in);
  ifstream eneFile(file.c_str(), ios::binary);

  if (!eneFile) {
    cerr << "The file can't be open." << endl;
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
	while (enVal <= 0.0) {
		enVal = bI.orbEne[ignoreInt];
		if (enVal < 0.0)
			ignoreInt++;
	};  

	bI.occOrb = ignoreInt;
	bI.virOrb = bI.rdim - bI.occOrb;

} 

void relGeometry::geomLecture(fstream& fileXYZ) {

//  string file = "../h2o.r4dscf.AtomUT";
  string file = "../h2.r4dscf.AtomUT";
  ifstream geomFile(file.c_str(), ios::in);

  if (!geomFile) {
    cerr << "The file can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  geomFile >> atmNum;
  geoRel.resize(atmNum);
  cout << "There are: " << atmNum << endl;
  for (int i=0; i<atmNum; i++)
    geomFile >> geoRel[i].atomType >> geoRel[i].atomicNum >> geoRel[i].nucCoor[0] >> geoRel[i].nucCoor[1] >> geoRel[i].nucCoor[2] >> geoRel[i].charge;

  geomFile.close();

  for (int i=0; i<atmNum; i++)
    cout << geoRel[i].atomType << " " <<  geoRel[i].atomicNum << " " << geoRel[i].nucCoor[0] << " " << geoRel[i].nucCoor[1] << " " << geoRel[i].nucCoor[2] << " " << geoRel[i].charge << endl;
  
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

//  string file = "../salidatesth2o";
  string file = "../salidatest";

  ifstream input;
 // input.open(file.c_str(), ios::binary);
  input.open(file.c_str(), ios::in);

  if (!input.is_open()) {
    cerr << "The file can't be open." << endl;
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

/*
  input.read((char*)&ignoreInt, 4);
  
	for (int j=0; j<fac; j++) {
		for (int k=0; k<fac; k++) { 
   		input.read((char*)&orb0, 8);
   		input.read((char*)&orb1, 8);
   		input.read((char*)&orb2, 8);
   		input.read((char*)&orb3, 8);
	 		bI.mo[0].moS[j].moRe[k] = orb0;
	 		bI.mo[1].moS[j].moRe[k] = orb1;
	 		bI.mo[0].moS[j].moIm[k] = orb2;
	 		bI.mo[1].moS[j].moIm[k] = orb3;
 	//		cout  << j << " " << k << " " << orb0 << " " << orb1 << "  " << orb2 << "  " << orb3 << endl;
			}
	}
	*/		
 input.close();

}


int main () {

  fstream xyz;
  fstream info;
  fstream basis;
  fstream coeffs;
	fstream ene;
  
  relBasis rB;
  rB.readBasisRel(info, basis, xyz, coeffs, ene);
//	rB.source();
  
  return 0;
}
