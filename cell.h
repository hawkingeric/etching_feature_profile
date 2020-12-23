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

        protected:
        private:
};




