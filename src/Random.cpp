#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Random.h"
using namespace std;

void RandomWalk::random() {

  string str;
  stringstream sstr;
  random_device r1;
  sstr << r1();
  sstr >> str;

  seed_seq seed1(str.begin(), str.end());
  g1.seed(seed1);
}

double RandomWalk::randomNumber(double value) {
  random_device rd;
  default_random_engine gen(rd());
  uniform_real_distribution<double>variant(0, value);

  return variant(gen);

}

double RandomWalk::uniform() {
	return uniformDistribution(g1);
}

double RandomWalk::normal(double mu, double sigma) {
  return sigma * normalDistribution(g1) + mu;
}

double RandomWalk::uniform (double low, double high) {
  return (high - low) * uniform() + low;
}

void RandomWalk::initialStep(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2) {

  double ran1, ran2, var, thetha;

  var = 1.0 - uniform();
  ran1 = sqrt(-2.0*log(var));
  thetha = 2.0*M_PI*uniform();
  var = 1.0 - uniform();
  ran2 = sqrt(-2.0*log(var));
  
  x1 += ran1*cos(thetha);
  y1 += ran1*sin(thetha);
  z1 += ran2;
	
  var = 1.0 - uniform();
  ran1 = sqrt(-2.0*log(var));
  thetha = 2.0*M_PI*uniform();
  var = 1.0 - uniform();
  ran2 = sqrt(-2.0*log(var));

  x2 += ran1*cos(thetha);
  y2 += ran1*sin(thetha);
  z2 += ran2;

}

/*
void RandomWalk::initialStep(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2) {

  double ran1, ran2, var, thetha1, thetha2;

  var = uniform() * 0.2;
  ran1 = sqrt(-0.5*log(var));
  var = uniform()*0.5;
  ran2 = sqrt(-0.5*log(var));
  thetha1 = uniform()*4.0*M_PI;
  thetha2 = uniform()*2.0*M_PI;

  x1 += ran1*cos(thetha1);
  y1 += ran1*sin(thetha1);
  z1 += ran2*cos(thetha2);

  var = uniform() * 0.2;
  ran1 = sqrt(-0.5*log(var));
  var = uniform()*0.5;
  ran2 = sqrt(-0.5*log(var));
  thetha1 = uniform()*4.0*M_PI;
  thetha2 = uniform()*2.0*M_PI;

  x2 += ran1*cos(thetha1);
  y2 += ran1*sin(thetha1);
  z2 += ran2*cos(thetha2);

}
*/

