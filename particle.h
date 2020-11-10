#include "etching.h"

double Maxwell_Boltzmann_pdf(double mass, double T, double v);
double Maxwell_Boltzmann_cdf(double mass, double T, double v);
double CalculateSpeedCutoff(double mass, double T);

class particle
{
        public:
                void setInitialPosition(double, double, double, double, double, double);
                double setInitialSpeed(double, double, double);
                double setInitialTheta_by_Gaussian(double);
                double setReemittedTheta(double);
                void setInitialType(std::vector<double>&);
                void ReflectedWithNewEnergy(double*, double*, double*, double*, double*, double*);
                void ReemittedWithNewDirection(double* , double, double, double);
                double dPos [3];
                double Vel [3];
                int iPos [3];
                double mass;
                double speed;
                double energy;
                double PropagationTimeInterval;
                double theta;
                double phi;
                int itagPos;
                int ParticleType;
                int iNumHitMask;
};
