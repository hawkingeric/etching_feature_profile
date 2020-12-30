#include "etching.h"
using namespace std;


class cell
{
        public:

                int iNumCell, iDimSize[3];
                double dDimLength[3], dDimLength_per_cell[3];
                int **iDim;
                int **iID_NBR;
                int *iStatus; // 0:vacuum 1:material 2: mask
                double *dNumCharge, *dNumMaterial, *dNumNeutral, *dNumMask, *dNumIon, **dNumSiClxg, **dNumSiClxs;
                double **ElectricForce;
                double dCoverage;
                double DeltaX, DeltaY, DeltaZ;


                void initial(int, int, int, double, double, double);
                void setNBR();
                void setStatus(int, int, int);
                void write_to_vtk(std::string, int, int, int, double, int*, double*, double*, double**, int, bool, bool, bool, bool, std::string, std::string, int);

                void Generation(double*, double*, double*,
                                    vector<double>, vector<double>, vector<double>,
                                    vector<double>, vector<double>, bool, bool, double, double, int*, double*, double*, int*, double*, double*, double* );
                        void GenerateParticleTypeMass(double*, int*, double*);
                        void GeneratePosition( double*, int*);
                        void GenerateTheta(double*);
                        void GeneratePhi(double*);
                        void GenerateRadicalSpeed(double*, double*, double*, double*);
                        void GenerateIonThetaEnergy(vector<double>&, vector<double>&,  vector<double>&, double*, double*);
                        void CalculateVelFromThetaPhi(double*, double*, double*, double*);
                        double Maxwell_Boltzmann_pdf(double, double, double);
                        double Maxwell_Boltzmann_cdf(double, double, double);
                        double CalculateSpeedCutoff(double, double);

                void Propagation(double*, int*, double*, int [][3], int*, int*, double*, double*, int*, double*, double*, double*);
                        void BoundaryMapping(int*, double*, int*);
                        void ParticleEdge(double*, double*, int [][3], int*);

                void SurfaceNormal(int [][3], int*,  int*, double*, double*, double*);
                        void SurfaceSites(int [][3] , int*, int*, double [][3], int* );

                void normReflectedVelocity(double*, double*, double*, double*, double*);
                void SurfaceAngle(double*, int*, double*, double*);

                void MaskReaction( int*, double*, double*, int*, double*);
                void SubstrateReaction( int*, int*, double*, double*, int*, double*);
                        void RadicalReactionProb(int*, int*, int*, double*, double* );
                        void IonReactionProb(int*, int*,  int*, double*, double*, double*, int*, double*, double*, double*, double*, double*);
                        void ClRadicalReaction(int*, double*, int*, int*);
                        void ClIonReaction(int*, double*, double*, double*, int*, int*, int*, int*);
                        void Cl2IonReaction(int*, double*, double*, double*, int*, int*, int*, int*);
                        void ArIonReaction(int*, double*, double*, double*, int*, int*, int*, int*);
                        void Redeposition( int*, int*, double*, int*, int*);


                void MaskReflection( double*, int*, double*, double*);
                void SubstrateReflection( double*, double*, double*, double*, int*, double*, double*, double*, double*);
                        //--Used in ProductReflection
                        void InelasticReflection(double*, double*, double*, double*, double*, double*);
                        void ElasticReflection(double*, double*, double*);
                        void DiffusiveReflection(double*, double*, double*, double*, double*, double*, double*);

        protected:
        private:
};




