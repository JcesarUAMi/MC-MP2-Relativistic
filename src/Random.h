#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <random>
#include <cstdio>
using namespace std;

class RandomWalk {
  public:
		void random();
    double randomNumber(double);
    double uniform(double, double);
    double uniform();
    double normal(double, double);
    void initialStep (double&, double&, double&, double&, double&, double&);
	private:
		mt19937 g1;
    normal_distribution<double> normalDistribution;
		uniform_real_distribution<double> uniformDistribution;

};

#endif


