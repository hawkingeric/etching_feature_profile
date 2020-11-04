#include "etching.h"

double Maxwell_Boltzmann_pdf(double mass, double T, double v);
double Maxwell_Boltzmann_cdf(double mass, double T, double v);
double calc_v_cut(double mass, double T);

class particle
{
        public:
                void setInitialPosition(double, double, double, double, double, double);
                double setInitialSpeed(double, double, double);
                double setInitialTheta_by_Gaussian(double);
                double setReemitTheta(double);
                void setInitialType(std::vector<double>&);
                void reflected_velocity_with_new_energy(double*, double*, double*);
                double dPos [3];
                double Vel [3];
                int iPos [3];
                double mass;
                double speed;
                double energy;
                double time_interval;
                double theta;
                double phi;
                int itagPos;
                int ParticleType;
                int iNumHitMask;
};
