#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include "General.h"
using namespace std;

void General::Data(fstream& basisFile, fstream& xyzFile, fstream& movecsFile) {

  readData.coefficientsFinal(basisFile, xyzFile, movecsFile);
  readData.nucleusBasis();
}

void General::OverlapMatrix() {
  
  int k, m, n, j, l, i, g, p;
  int ang1[3], ang2[3]; 
  double Ax, Bx, Px, Ay, By, Py, Az, Bz, Pz;;
  double PAx, PBx, PAy, PBy, PAz, PBz, AB;
  double expo1, expo2, mu, angMom1, angMom2, expoTot, funcType1, funcType2;
  double traslapeTot, Kab, traslapeX, traslapeY, traslapeZ;
  int num = 0;
  for (m=0; m<readData.lastBasis.size(); m++)
    num += readData.lastBasis[m].basisExponents.size();
 
  overlapElements.resize(num*num); 
  
  m = 0;
  for (k=0; k<readData.lastBasis.size(); k++) {
    i = 0;
    for (n=0; n<readData.lastBasis[k].basisExponents.size(); n++) {
      angMom1 = readData.lastBasis[k].angularMoment[n];
      funcType1 = readData.lastBasis[k].functionType[n];
      expo1 = readData.lastBasis[k].basisExponents[n];
      basis.angularMoment(ang1, angMom1, funcType1);
      for (j=0; j<readData.lastBasis.size(); j++) {
        g = 0; 
        for (l=0; l<readData.lastBasis[j].basisExponents.size(); l++) {
          angMom2 = readData.lastBasis[j].angularMoment[l];
          funcType2 = readData.lastBasis[j].functionType[l];
          expo2 = readData.lastBasis[j].basisExponents[l];
          basis.angularMoment(ang2, angMom2, funcType2);

          expoTot = expo1 + expo2;
          mu = (expo1 * expo2) / expoTot;

          Ax = readData.lastBasis[k].nucleusCoord[i];
          Ay = readData.lastBasis[k].nucleusCoord[i+1];
          Az = readData.lastBasis[k].nucleusCoord[i+2];

          Bx = readData.lastBasis[j].nucleusCoord[g];
          By = readData.lastBasis[j].nucleusCoord[g+1];
          Bz = readData.lastBasis[j].nucleusCoord[g+2];

          Px = (expo1 * Ax + expo2 * Bx) / expoTot; 
          Py = (expo1 * Ay + expo2 * By) / expoTot; 
          Pz = (expo1 * Az + expo2 * Bz) / expoTot; 
         
          AB = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);

          Kab = -mu * AB;

          if (abs(Kab) >= 45) 
            traslapeTot = 0.0;
          else {
            Kab = exp(Kab);

            PAx = Px - Ax;
            PAy = Py - Ay;
            PAz = Pz - Az;

            PBx = Px - Bx;
            PBy = Py - By;
            PBz = Pz - Bz;

            traslapeX = TraslapeNumericalIntegral(ang1[0], ang2[0], PAx, PBx, expoTot);
            traslapeY = TraslapeNumericalIntegral(ang1[1], ang2[1], PAy, PBy, expoTot);
            traslapeZ = TraslapeNumericalIntegral(ang1[2], ang2[2], PAz, PBz, expoTot);
            traslapeTot = Kab * traslapeX * traslapeY * traslapeZ;
/*
          cout << expo1 << " " << Ax << " " << Ay << " " << Az << " angular moments " << ang1[0] << " " << ang1[1] << " " << ang1[2] << endl;
          cout << expo2 << " " << Bx << " " << By << " " << Bz << " angular moments " << ang2[0] << " " << ang2[1] << " " << ang2[2] << endl;
          cout << expoTot << " " << Px << " " << Py << " " << Pz << endl;
          cout << traslapeTot << "\n\n";*/
          }
        overlapElements[m] = traslapeTot;
        m++;
        }
      }
    }  
  }
}

double General::TraslapeNumericalIntegral(int la, int lb, double Xpa, double Xpb, double q) {
    int i;
    double integral = 0.0, sq;
    double x, xk[12], wk[12];
    sq = sqrt(q);

  xk[0]  = -3.889724897869781919272;     xk[1]  = -3.020637025120889771711;
  xk[2]  = -2.279507080501059900188;     xk[3]  = -1.59768263515260479671;
  xk[4]  = -0.9477883912401637437046;    xk[5]  = -0.314240376254359111277;
  xk[6]  =  3.889724897869781919272;     xk[7]  =  3.020637025120889771711;
  xk[8]  =  2.279507080501059900188;     xk[9]  =  1.59768263515260479671;
  xk[10] = 0.9477883912401637437046;     xk[11] =  0.314240376254359111277;
  wk[0] = 2.65855168435630160602E-7;     wk[1]  = 8.5736870435878586546E-5;
  wk[2] = 0.00390539058462906185999;     wk[3]  = 0.05160798561588392999187;
  wk[4] = 0.2604923102641611292334;      wk[5]  = 0.5701352362624795783471;
  wk[6] = 2.65855168435630160602E-7;     wk[7]  = 8.5736870435878586546E-5;
  wk[8] = 0.00390539058462906185999;     wk[9]  = 0.05160798561588392999187;
  wk[10] = 0.2604923102641611292334;     wk[11] = 0.5701352362624795783471;

  if (la == 0 && lb == 0) 
    return sqrt(M_PI) / sq;
  else {
    for (i = 0; i < 12; i++) {
      x = xk[i]/sq;
      integral += wk[i] * pow(x + Xpa, la) * pow(x + Xpb, lb);
    }
  return integral/sq;
  }

}
