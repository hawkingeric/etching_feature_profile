#include "etching.h"
using namespace std;


class cell
{
        public:
                void initial(int, int, int, double, double, double);
                void setNBR();
                void setStatus(int, int, int);
                int iNumCell, iDimSize[3];
                double dDimLength[3], dDimLength_per_cell[3];
                int **iDim;
                int **iID_NBR;
                int *iStatus; // 0:vacuum 1:material 2: mask
                double *dNumCharge, *dNumMaterial, *dNumNeutral, *dNumMask, *dNumIon, **dNumSiClxg, **dNumSiClxs;
                double **ElectricForce;
                double dCoverage;
                double DeltaX, DeltaY, DeltaZ;
                void Generation(int, double, double*, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&,
                                                    bool, bool, double, double, int*, double*, double**, double**, double*);
                void Propagation(int, double, int, double, double, int, int*, double*, double**, double**, double* );
                void BoundaryMapping(int*, double*, int*);
                void CalcParticleEdge( double*, int*, int*,  double [][3], double, double, double, double [][3], int [][3]  );
                void InelasticReflection(double*, double*, double*, double, double*);
                void ElasticReflection(double*, double, double*);
                void DiffusiveReflection(double, double*, double, double*, double, double*);
                void CalcSurfaceNormal(int [][3], int, int, int*, double*, double*, double*, double*, double* );
                void SurfaceReaction( int , int, double, double, int*, double*, double);
                void ClRadicalReaction(int, int*, int*);
                void ClIonReaction(int, double, double, int*, int*, int*);
                void Cl2IonReaction(int, double, double, int*, int*, int*);
                void ArIonReaction(int, double, double, int*, int*, int*);
                void Redeposition( int, int, int*, int*);
        protected:
        private:
};




