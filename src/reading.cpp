#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
using namespace std;

int main () {

  int N, n, Nbs, Nb; 
  double Nbi;
  cout << "Introduce the number of steps on total: ";
  cin >> N;
  cout << "Introduce the size of the block: ";
  cin >> Nbs;
  Nb = N/Nbs;
  Nbi = 1.0/Nb;
  vector<double> result;
  vector<double> emp2Var;
  result.resize(Nb);
  emp2Var.resize(Nb);
  double finalVar;
  double finalRes;
  int j;

  string file = "testH210";
//  string file = "testH2Om12";
//  string file = "testCuHm12";
//  string file = "testCu2m12";
//  string file = "testAgH";
//  string file = "testAg2m8";
//  string file = "testAu2m12";
//  string file = "testAuHm12";

  ifstream data(file.c_str(), ios::in);

  if (!data) {
    cerr << "The MC file can't be open." << endl;
    exit(EXIT_FAILURE);
  }

  finalRes = 0.0;
  finalVar = 0.0;
  for (int i=0; i<Nb; i++) {
    data >> result[i] >> emp2Var[i];
    finalRes += result[i]; 
    finalVar += emp2Var[i];
    j = i+1;
//   cout << i << "      " << /*(-finalRes/j - (0.321182739))*1000  << "    " <<  sqrt(Nbi*finalVar/double(j*(j-1.0)))*1000 << endl;
  }

  data.close();

  double other;
  other = sqrt(Nbi*finalVar/double(Nb*(Nb-1.0)));
  cout << "Aprox:  " << finalRes/Nb << " with:   "  <<  other << endl;


 return 0; 

}
