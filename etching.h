#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <random>
#include <chrono>
#include <math.h>
#include <vector>

#define PI acos(-1)
#define Joule_to_eV 6.24150913E18
#define BoltzmannConstant_SI   1.38064852E-23  // unit: Joule*K^-1
#define BoltzmannConstant    8.617333262145E-5 // unit: eV*K^-1

void write_to_vtk(std::string, int, int, int, double, int*, double*, double*, double**, int, bool, bool, bool, bool, std::string, std::string, int);
std::vector<double> gauss(std::vector<std::vector<double>> A);
std::vector<double> read_data(std::string);


void Generation(int, double, double, double, double, double, double, double, double*,
                                    std::vector<double>&, std::vector<double>&, std::vector<double>&,
                                    std::vector<double>&, std::vector<double>&,
                                    bool, bool,
                                    double, double,
                                    int*, double*, double**, double*, double**);

void Propagation(int, double, int, double,
                                      int*, double*, double**, double*, double**,
                                      double, double, double, int, int, int, double, double, double,
                                      double [][3], double, double, double,
                                      int*, int**, double*, double**, double**,
                                      int [][3] , int);
void ReflectedWithNewEnergy(double*, double*, int, double*, double*, double*, double*);
void ReemittedWithNewDirection(double*, double, double, double, double*);
void surface_normal(int [][3], int, int, int*, double*, int, int, int, int*, int**, int, int, int, double*, double*, double*, double* );
void surface_reaction( int , int, double, double, int*, double*, double*, int*, double*, double**);
void ClRadicalReaction(int, int*, double**, double*, int*, int*);
void ClIonReaction(int, int*, double**, double*, double, double, int*, int*, int*, int*);
void Cl2IonReaction(int, int*, double**, double*,  double, double, int*, int*, int*, int*);
void ArIonReaction(int, int*, double**, double*, double, double, int*, int*, int*, int*);
void redeposition( int, int, int*, double**, double*, int*, int*);

