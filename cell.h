#include "etching.h"
using namespace std;
class cell
{
        public:
                void initial(int, int, int, double, double, double);
                void setNBR();
                void setStatus(int, int, int);
                void setStatusAll(int, int);
                void setStatusPlane(int, int, int, int);
                void setStatusLine(int *, int *, int, int, int);
                void setInterFace(int);
                void setInterFace2(int);
                void setInterFace4Vacuum(int i);
                void TransNeutral4Vacuum(int i);
                void printID(int);
                void printNBR(int);
                void Slope(int, double *);
                void Slope2(int, double *);
                void Slope3(int, double *);
                void surface_normal(int [][3] , int, int, int*, double*, double*, double*, double*, double*);
                void ClRadicalReaction(vector<double>&, int, int, int*, int*);
                void redeposition(vector<double>&, int, int, int*, int*);
                void ArIonReaction(vector<double>&, double*, vector<double>&, vector<int>&, vector<double>&, vector<double>&, int,  int, double, double, int*, int*, int*, int*);
                void ClIonReaction(vector<double>&, double*, vector<double>&, vector<int>&, vector<double>&, vector<double>&,  int, int, double,  double, int*, int*, int*, int*);
                void Cl2IonReaction(vector<double>&, double*, vector<double>&, vector<int>&, vector<double>&, vector<double>&, int, int, double, double, int*, int*, int*, int*);
                void IonMaskReaction(vector<double>&, int, double,  double , int* );



                int iNumCell, iDimSize[3];
                double dDimLength[3], dDimLength_per_cell[3];
                int **iDim;
                int **iID_NBR;
                int *iStatus; // 0:vacuum 1:material 2: mask
                double *dNumCharge, *dNumMaterial, *dNumNeutral, *dNumMask, *dNumIon, **dNumSiClxg, **dNumSiClxs;
                double **grad_potential;
                double dCoverage;
                double DeltaX, DeltaY, DeltaZ;


        protected:
        private:
};




