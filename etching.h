#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <random>
#include <chrono>
#include <vector>

#define PI acos(-1)
#define Joule_to_eV 6.24150913E18
#define BoltzmannConstant_SI   1.38064852E-23  // unit: Joule*K^-1
#define BoltzmannConstant    8.617333262145E-5 // unit: eV*K^-1

extern double MassChlorine;
extern double MassArgon;
extern double MassElectron;
extern double MassSilicon;
extern std::ofstream debugfile;
extern int iVacuumStat, iSubstrateStat, iMaskStat, iInterS_VStat, iInterM_VStat;
extern const int iClRadicalType, iArIonType, iClIonType, iCl2IonType, iSigType, iSiClgType, iSiCl2gType, iSiCl3gType, iSiCl4gType ;
extern int X_dir, Y_dir, Z_dir;
extern int P_sputtering, C_sputtering;


void write_to_vtk(std::string, int, int, int, double, int*, double*, double**, int, bool, bool, bool, bool, std::string, std::string, int);
std::vector<double> gauss(std::vector<std::vector<double>> A);
std::vector<double> read_data(std::string);




