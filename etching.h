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

extern int X_dir, Y_dir, Z_dir;
extern int iClRadicalType, iSigType, iSiClgType, iSiCl2gType, iSiCl3gType, iClIonType, iCl2IonType, iArIonType;
extern int iVacuumStat, iSubstrateStat, iMaskStat;
extern double MassChlorine, MassArgon, MassSilicon, MassssElectron;
extern int P_sputtering, C_sputtering;

std::vector<double> gauss(std::vector<std::vector<double>> A);
std::vector<double> read_data(std::string);




