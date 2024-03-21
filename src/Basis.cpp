#include <iostream>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
using namespace std;

#include "Basis.h"

void Basis::readBasisFile(fstream& readFromFile, int nDifAtoms) {
  
//  string file = "../cc-pvdz.1_Be.gbs";
//  string file = "6-311++g_st.gbs";
//  string file = "../6-31g_st.gbs";
//  string file = "../6-31g.1.gbs";
//  string file = "sto-3g.1.gbs";
//  string file = "cc-pwcvtz.0.gbs";
    string file = "../cc-pvdz.1.gbs";
//    string file = "3-21g.1.gbs";

//  if (type == "Spherical") 
    SoC = 0; // spherical by default
//  else {  // for cartesian symmetry
//    SoC = 1;
//  }

  ifstream basis(file.c_str(), ios::in);
  
  if (!basis) {
    cerr << "The file can't be open." << endl;
    exit(EXIT_FAILURE);
  }
  string atom;
  string primitiveType;
  string lines;
  string lineNumbers;
  int ang;
  double exponent;
  double coeff;
  double *coeffSP;
  double *exponentSP;
  int symmetry;
  int primitiveNumber;
  basisPrimitiveFunctions.resize(nDifAtoms);
  int t;
  size_t position;
  size_t result;
  do {
    getline(basis, lines);
    t = lines.size();
  } while (t >= 12);
  getline(basis, lines);
  
  for (t=0; t<nDifAtoms; t++) {
    result = 10;
    getline(basis, lines);
    istringstream lineCheck(lines);
    lineCheck >> atom;
    Molecule::getAtomicNumber(atom);    
    basisPrimitiveFunctions[t].atomType = Molecule::atomicNumber;
  while (result > 5) {
    getline(basis, lineNumbers);
    result = lineNumbers.size();
    if (result > 4) { 
      istringstream lineCheck(lineNumbers);
      lineCheck >> primitiveType >> primitiveNumber;
      ang = getSymmetry(primitiveType, SoC);
      if (primitiveType == "SP") {
       coeffSP = new double[primitiveNumber]; 
       exponentSP = new double[primitiveNumber]; 
       basisPrimitiveFunctions[t].shellSize.push_back(primitiveNumber);
        for (int i = 0; i<primitiveNumber; i++) {
          getline(basis, lineNumbers);
          position = lineNumbers.find("D");
          while(position != string::npos) {
            lineNumbers.replace(position, 1, "E");
            position = lineNumbers.find("D", position+1);
              }
          istringstream lineCheck(lineNumbers);
          lineCheck >> exponentSP[i] >> coeff >> coeffSP[i];
          basisPrimitiveFunctions[t].primitiveExponents.push_back(exponentSP[i]);
          basisPrimitiveFunctions[t].primitiveCoefficients.push_back(coeff);
          basisPrimitiveFunctions[t].totalAngMom.push_back(0);
          basisPrimitiveFunctions[t].primitiveTypes.push_back(0);
        }
        for (int j=0; j<3; j++) {
          basisPrimitiveFunctions[t].shellSize.push_back(primitiveNumber);
          for (int in=0; in<primitiveNumber; in++) {
            basisPrimitiveFunctions[t].totalAngMom.push_back(1);
            basisPrimitiveFunctions[t].primitiveExponents.push_back(exponentSP[in]);
            basisPrimitiveFunctions[t].primitiveCoefficients.push_back(coeffSP[in]);
            if (SoC == 0)
              basisPrimitiveFunctions[t].primitiveTypes.push_back(j);
            else 
              basisPrimitiveFunctions[t].primitiveTypes.push_back(j);
            
          }
        }
      delete[] coeffSP;
      delete[] exponentSP;
      } else {
        if (SoC == 0)
          symmetry = 2*ang + 1;
        else
          symmetry = (ang+2)*(ang+1)/2;

        coeffSP = new double[primitiveNumber]; 
        exponentSP = new double[primitiveNumber]; 
        for (int i = 0; i<primitiveNumber; i++) {
          getline(basis, lineNumbers);
          position = lineNumbers.find("D");
          while(position != string::npos) {
            lineNumbers.replace(position, 1, "E");
            position = lineNumbers.find("D", position+1);
              }
            istringstream lineCheck(lineNumbers);
            lineCheck >> exponentSP[i] >> coeffSP[i];
          }
        for (int j=0; j<symmetry; j++) {
          basisPrimitiveFunctions[t].shellSize.push_back(primitiveNumber);
          for (int i = 0; i<primitiveNumber; i++) {
            basisPrimitiveFunctions[t].totalAngMom.push_back(ang);
            basisPrimitiveFunctions[t].primitiveExponents.push_back(exponentSP[i]);
            basisPrimitiveFunctions[t].primitiveCoefficients.push_back(coeffSP[i]);
            if (SoC == 0)
              basisPrimitiveFunctions[t].primitiveTypes.push_back(j);
            else
              basisPrimitiveFunctions[t].primitiveTypes.push_back(j); 
              }
            }
        delete[] coeffSP;
        delete[] exponentSP;
          }
        }
      }
    }

  
  basis.close();

  for (int i=0; i<nDifAtoms; i++) {
    basisPrimitiveFunctions[i].normConst.resize(basisPrimitiveFunctions[i].primitiveExponents.size());
    for (int u=0; u<basisPrimitiveFunctions[i].normConst.size(); u++) 
      basisPrimitiveFunctions[i].normConst[u] = 0.0;
  }

  normalize(nDifAtoms);


  for(int k=0; k<nDifAtoms; k++) {
    cout << " another atom dude" << endl;
    for(int j=0; j<basisPrimitiveFunctions[k].shellSize.size(); j++) {
      cout << j << " " << basisPrimitiveFunctions[k].totalAngMom[j] << " size:   " << basisPrimitiveFunctions[k].shellSize[j] << endl;
    }
  }
 
  for (int k=0; k<nDifAtoms; k++) {
    cout << basisPrimitiveFunctions[k].atomType << " challe " << endl;
    cout << "Exponent " << '\t' << "Coefficient" << '\t' << "Type" << endl;
    for (int i=0; i<basisPrimitiveFunctions[k].primitiveExponents.size(); ++i) {
      cout << i << "             " << basisPrimitiveFunctions[k].primitiveExponents[i] << '\t' << basisPrimitiveFunctions[k].primitiveCoefficients[i] << '\t' << "\t " << basisPrimitiveFunctions[k].totalAngMom[i] << " \t" << basisPrimitiveFunctions[k].primitiveTypes[i] << endl;
      }
  }

}

int Basis::getSymmetry (string gT, int SoC) {

  if (gT == "S") 
      return 0;
   else if (gT == "P") 
      return 1;
   else if (gT == "D") 
      return 2;
   else if (gT == "F") 
      return 3;
   else if (gT == "G") 
      return 4;
   else if (gT == "H") 
      return 5;
   else if (gT == "I") 
      return 6;
}

void Basis::normalize(int nDifAtoms) {

  double coeff;
  double expo;
  int angMom;
  int l;
  int type;
  int vang[3];

  if (SoC == 0) {
    for (int k=0; k<nDifAtoms; k++) {
      l = 0;
      for (int i=0; i<basisPrimitiveFunctions[k].shellSize.size(); i++) {
        for (int j=0; j<basisPrimitiveFunctions[k].shellSize[i]; j++) {
          expo = basisPrimitiveFunctions[k].primitiveExponents[l];
          angMom = basisPrimitiveFunctions[k].primitiveTypes[l];
          type = basisPrimitiveFunctions[k].totalAngMom[l];
          basisPrimitiveFunctions[k].normConst[l] = normalizeSpherical(expo, type);
          basisPrimitiveFunctions[k].primitiveCoefficients[l] *= basisPrimitiveFunctions[k].normConst[l];
          l++;
        }
      }
    }
  } else {
    for (int k=0; k<nDifAtoms; k++) {
      for (int j=0; j<basisPrimitiveFunctions[k].primitiveExponents.size(); j++) {
        angMom = basisPrimitiveFunctions[k].totalAngMom[j];
        type = basisPrimitiveFunctions[k].primitiveTypes[j];
        angularMoment(vang, angMom, type);
        expo = basisPrimitiveFunctions[k].primitiveExponents[j];
        basisPrimitiveFunctions[k].normConst[j] = normalizeCartesian(expo, vang);
        basisPrimitiveFunctions[k].primitiveCoefficients[j] *= basisPrimitiveFunctions[k].normConst[j]; 
      }
    }
  }
}


void Basis::angularMoment (int *vang, int ang, int type) {

  for (int i=0; i<3; i++)
    vang[i] = 0;

  if (ang == 0) {
    vang[0] = 0;
    vang[1] = 0;
    vang[2] = 0;
  } else if  (ang == 1) {
    switch (type) {
      case 0 : vang[0]=1;  
               break;
      case 1 : vang[1]=1;  
               break;
      case 2 : vang[2]=1;  
               break;
      default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 2) {
    switch (type) {
      case 0 : vang[0]=2;  
               break;
      case 1 : vang[0]=1;
               vang[1]=1;
               break;
      case 2 : vang[1]=1;
               vang[2]=1;  
               break;
      case 3 : vang[1]=2;  
               break;
      case 4 : vang[0]=1;
               vang[2]=1;  
               break;
      case 5 : vang[2]=2;  
               break;
      default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 3) {
    switch (type) {
      case 0: vang[0]=3;  
              break;
      case 1: vang[0]=2;
              vang[1]=1;  
              break;
      case 2: vang[0]=2;
              vang[2]=1;  
              break;
      case 3: vang[0]=1;
              vang[1]=2;  
              break;
      case 4: vang[0]=1;
              vang[1]=1;
              vang[2]=1;  
              break;
      case 5: vang[0]=1;
              vang[2]=2;  
              break;
      case 6: vang[1]=3;
              break;
      case 7: vang[1]=2;
              vang[2]=1;  
              break;
      case 8: vang[1]=1;
              vang[2]=2;  
              break;
      case 9: vang[2]=3;
              break;
      default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 4) {
    switch (type) {
      case 0: vang[0]=4;  
               break;
      case 1: vang[0]=3;
              vang[1]=1;  
              break;
      case 2: vang[0]=3;
              vang[2]=1;  
              break;
      case 3: vang[0]=2;
              vang[1]=2;  
              break;
      case 4: vang[0]=2;
              vang[1]=1;
              vang[2]=1;  
              break;
      case 5: vang[0]=2;
              vang[2]=2;  
              break;
      case 6: vang[0]=1;
              vang[1]=3;  
              break;
      case 7: vang[0]=1;
              vang[1]=2;
              vang[2]=1;  
              break;
      case 8: vang[0]=1;
              vang[1]=1;
              vang[2]=2;  
              break;
      case 9: vang[0]=1;
              vang[2]=3;  
              break;
      case 10: vang[1]=4;
               break;
      case 11: vang[1]=3;
               vang[2]=1;  
               break;
      case 12: vang[1]=2;
               vang[2]=2;  
               break;
      case 13: vang[1]=1;
               vang[2]=3;  
               break;
      case 14: vang[2]=4;
               break;
    default: cout << "Check the type of primitive." << endl;
    }
  } else if (ang == 5) {
    switch (type) {
      case 0: vang[1]=5;  
              break;
      case 1: vang[1]=1;
              vang[2]=4;  
              break;
      case 2: vang[1]=2;
              vang[2]=3;  
              break;
      case 3: vang[1]=3;
              vang[2]=2;  
              break;
      case 4: vang[1]=4;
              vang[2]=1;  
              break;
      case 5: vang[2]=5;  
              break;
      case 6: vang[0]=1;
              vang[2]=4;  
              break;
      case 7: vang[0]=1;
              vang[1]=1;
              vang[2]=3;  
              break;
      case 8: vang[0]=1;
              vang[1]=2;
              vang[2]=2;  
              break;
      case 9: vang[0]=1;
              vang[1]=3;
              vang[2]=1;  
              break;
      case 10: vang[0]=1;
               vang[1]=4;  
               break;
      case 11: vang[0]=2;
               vang[2]=3;  
               break;
      case 12: vang[0]=2;
               vang[1]=1;
               vang[2]=2;  
               break;
      case 13: vang[0]=2;
               vang[1]=2;
               vang[2]=1;  
               break;
      case 14: vang[0]=2;
               vang[1]=3;  
               break;
      case 15: vang[0]=3;
               vang[2]=2;  
               break;
      case 16: vang[0]=3;
               vang[1]=1;
               vang[2]=1;  
               break;
      case 17: vang[0]=3;
               vang[1]=2;  
               break;
      case 18: vang[0]=4;
               vang[2]=1;  
               break;
      case 19: vang[0]=4;
               vang[1]=1;  
               break;
      case 20: vang[0]=5;  
               break;
      default: cout << "Check the type of primitive." << endl;
    }
  }
}


double Basis::normalizeSpherical (double exponent, int ang) {
  double result;
  result = pow(2.0*exponent/M_PI,0.75);
  result *= sqrt(pow(4*exponent, ang) / factorialDouble(2*ang-1));

/*  result1 = 2.0*pow(2.0*exponent, ang + 1.5);
  factor = tgamma(ang + 1.5);
  result1 = sqrt(result1/factor);
*/
  cout << "ya pues: " << exponent << "   " << ang << "   " << result << endl;


  return result;
}

double Basis::normalizeCartesian (double exponent, int *ang) {

  double result;
  double finalResult;
  double pi = 2*exponent/M_PI;
  finalResult = 1.0;
  for(int i=0; i<3; i++) {
    result = pow(pi, 0.25);
    if (ang[i] != 0) { 
      result *= pow(4*exponent, 0.5*ang[i]);
      result /= sqrt(factorialDouble(2*ang[i] - 1));
    }
    finalResult *= result;
  }

  return finalResult;
}


double Basis::factorialDouble (int k) {
  if (k <= 0)
    return 1.0;
  else
    return k*factorialDouble(k-2);
/*
  int value;
  int result;
  result = 1;
  if (n > 0) {
    value = n/2 - 1;
    for (int i=0; i<=value; i++)
      result *= (n - 2*i);
  }
  
 return result; 
*/
 }
