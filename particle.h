#include "etching.h"

double Maxwell_Boltzmann_pdf(double mass, double T, double v);
double Maxwell_Boltzmann_cdf(double mass, double T, double v);
double CalculateSpeedCutoff(double mass, double T);

class particle
{
        public:
                void setInitialPosition(double, double, double, double, double, double);
                void setInitialType(double*);
                double setInitialSpeed(double, double, double);

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
