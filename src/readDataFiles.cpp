#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "readDataFiles.h"
using namespace std;

void readDataInput::coefficientsFinal(fstream& basisFile, fstream& xyzFile, fstream& movecsFile) {

//  molec.readCoordinates(xyzFile);
  molec.getNumberOfDifAtoms();
  lastBasis.resize(molec.nNuc);
  basis.readBasisFile(basisFile, molec.nDifAtoms);
  coeff.readMovecsFile(movecsFile);
}

void readDataInput::nucleusBasis() {
  
  int atom;
  int j;
  numberOfMolOrb = coeff.c.nw_nbf;

  for (size_t i=0; i<molec.nNuc; i++) {
    atom = molec.geom[i].nAtomic;
    j = 0;
    while (j<molec.nDifAtoms) {
      if (atom == basis.basisPrimitiveFunctions[j].atomType) {
        lastBasis[i].atomType = atom;
        lastBasis[i].shellSize.resize(basis.basisPrimitiveFunctions[j].shellSize.size());
        lastBasis[i].functionType.resize(basis.basisPrimitiveFunctions[j].primitiveTypes.size());
        lastBasis[i].angularMoment.resize(basis.basisPrimitiveFunctions[j].totalAngMom.size());
        lastBasis[i].basisCoefficients.resize(basis.basisPrimitiveFunctions[j].primitiveCoefficients.size());
        lastBasis[i].basisExponents.resize(basis.basisPrimitiveFunctions[j].primitiveExponents.size());
        copy(begin(molec.geom[i].nCoor), end(molec.geom[i].nCoor), begin(lastBasis[i].nucleusCoord));  
        copy(begin(basis.basisPrimitiveFunctions[j].shellSize), end(basis.basisPrimitiveFunctions[j].shellSize), begin(lastBasis[i].shellSize));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveExponents), end(basis.basisPrimitiveFunctions[j].primitiveExponents), begin(lastBasis[i].basisExponents));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveCoefficients), end(basis.basisPrimitiveFunctions[j].primitiveCoefficients), begin(lastBasis[i].basisCoefficients));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveTypes), end(basis.basisPrimitiveFunctions[j].primitiveTypes), begin(lastBasis[i].functionType));
        copy(begin(basis.basisPrimitiveFunctions[j].totalAngMom), end(basis.basisPrimitiveFunctions[j].totalAngMom), begin(lastBasis[i].angularMoment));
        j = molec.nDifAtoms;
      } else
        j++;
    }
  }

  int cont, h;
  cont = 0;
  numberOfPrimFunc = 0;
  
  moisCoeffs.resize(numberOfMolOrb);
  for (int j=0; j<numberOfMolOrb; j++) {
    moisCoeffs[j].resize(numberOfMolOrb);
    for (int i=0; i<numberOfMolOrb; i++) { 
      moisCoeffs[j][i] = coeff.molecularOrbitalData[j].orbitalCoeff[i];
    }
  }


  for (int i=0; i<lastBasis.size(); i++) {
    numberOfPrimFunc += lastBasis[i].basisExponents.size();  
    lastBasis[i].lastBasisCoeff.resize(coeff.c.nw_nbf); 
    for(int u=0; u<lastBasis[i].lastBasisCoeff.size(); u++) {
      lastBasis[i].lastBasisCoeff[u].bC.resize(lastBasis[i].shellSize.size());
      h = 0;
      for (int k=0; k<lastBasis[i].shellSize.size(); k++) {
        lastBasis[i].lastBasisCoeff[u].bC[k].resize(lastBasis[i].shellSize[k]);
        for (int l=0; l<lastBasis[i].shellSize[k]; l++) {
          lastBasis[i].lastBasisCoeff[u].bC[k][l] = coeff.molecularOrbitalData[u].orbitalCoeff[cont+k] * lastBasis[i].basisCoefficients[h];
           h++;
        }
      }
    }
   cont += lastBasis[i].shellSize.size();
  }

}

void readDataInput::densityMatrix() {
  
  int cont, cont1;

  getNumOccOrb();
  dEi.resize(numberOfOccOrb);
  for (int i=0; i<numberOfOccOrb; i++)
    for (int k=0; k<numberOfPrimFunc; k++) {
      dEi[i].dE.resize(numberOfPrimFunc);
      for (int j=0; j<numberOfPrimFunc; j++) {
        dEi[i].dE[k].resize(numberOfPrimFunc);
        dEi[i].dE[k][j] = 0.0;
    }
  }
 
  cont1 = 0; 
  for (int j=0; j<lastBasis.size(); j++) 
    for (int k=0; k<lastBasis[j].shellSize.size(); k++) 
      for (int l=0; l<lastBasis[j].shellSize[k]; l++) {
        cont = 0;
        for (int h=0; h<lastBasis.size(); h++) 
          for (int m=0; m<lastBasis[h].shellSize.size(); m++)
            for (int n=0; n<lastBasis[h].shellSize[m]; n++) { 
              for (int i=0; i<numberOfOccOrb; i++)  
              dEi[i].dE[cont][cont1] = 2.0*lastBasis[h].lastBasisCoeff[i].bC[m][n] * lastBasis[j].lastBasisCoeff[i].bC[k][l];    
              
                cont++;
              }
          cont1++;
        }
}  

void readDataInput::getNumOccOrb () {

  int occ;
  double energy;
  numberOfOccOrb = 0;
  numberOfVirOrb = 0; 
  for(int m=0; m<numberOfMolOrb; m++) {
    occ = coeff.molecularOrbitalData[m].orbitalOccupancy;
    energy = coeff.molecularOrbitalData[m].orbitalEnergie;
    if (occ != 0) {
      numberOfOccOrb++;
      Occ.push_back(energy);
    } else {
      numberOfVirOrb++;
      Vir.push_back(energy);
    }
  }
  numberOfPrimFunc = 0;
  for (int i=0; i<lastBasis.size(); i++) 
    numberOfPrimFunc += lastBasis[i].basisExponents.size();

  moi.resize(numberOfPrimFunc);

}

double readDataInput::generateMoiSpherical(double x, double y, double z) {
  
  double rr;
  double distx, disty, distz;
  double fac;
  double expo, prim, mo, den, coef;
  int ang[3];
  int angMom;
  int funcType;
  int k;
  k = 0;
  for (int i=0; i<lastBasis.size(); i++)
    for (int j=0; j<lastBasis[i].basisExponents.size(); j++) {
      distx = x - lastBasis[i].nucleusCoord[0];
      disty = y - lastBasis[i].nucleusCoord[1];
      distz = z - lastBasis[i].nucleusCoord[2];
      rr = distx*distx + disty*disty + distz*distz;
      expo = exp(-lastBasis[i].basisExponents[j]*rr);
      coef = lastBasis[i].basisCoefficients[k];
      angMom = lastBasis[i].angularMoment[j];
      funcType = lastBasis[i].functionType[j];
      fac = SphericalEvaluation(distx, disty, distz, rr, angMom, funcType);
      moi[k] = fac*expo;
      k++;
    }


  den = 0.0;
  for (int i=0; i<numberOfOccOrb; i++) {
    mo = 0.0;
    k = 0;
    for (int h=0; h<lastBasis.size(); h++)
      for (int m=0; m<lastBasis[h].shellSize.size(); m++)
        for (int n=0; n<lastBasis[h].shellSize[m]; n++) {
          prim = moi[k]*moi[k]*lastBasis[h].lastBasisCoeff[i].bC[m][n]* lastBasis[h].lastBasisCoeff[i].bC[m][n];
          mo += prim;
          den += prim;
          k++;
        }
      }

  return den;

}

double readDataInput::SphericalEvaluation(double x, double y, double z, double r, int ang, int type) {

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
      case 3 : result = sc[0]*x*z;
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

