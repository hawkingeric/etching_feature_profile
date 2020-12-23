#include "etching.h"
#include "cell.h"
#include "jacobi_eigenvalue.h"
#include "rand.h"

using namespace std;
#pragma acc routine seq
void cell::initial(int ix, int iy, int iz, double Lx, double Ly, double Lz)
{
        int X_dir = 0;
        int Y_dir = 1;
        int Z_dir = 2;

        cell::iDimSize[0] = ix;
        cell::iDimSize[1] = iy;
        cell::iDimSize[2] = iz;
        cell::dDimLength[0] = Lx;
        cell::dDimLength[1] = Ly;
        cell::dDimLength[2] = Lz;
        cell::dDimLength_per_cell[0]= Lx/ix;
        cell::dDimLength_per_cell[1] = Ly/iy;
        cell::dDimLength_per_cell[2] = Lz/iz;

        int i;
        cell::iNumCell       = cell::iDimSize[X_dir] * cell::iDimSize[Y_dir] * cell::iDimSize[Z_dir]; // ix * iy * iz
        cell::iStatus        = new int [cell::iNumCell];    // 1-D array int type with size of iNumCell
        cell::dNumCharge     = new double [cell::iNumCell]; // 1-D array double type with size of iNumCell
        cell::dNumMaterial   = new double [cell::iNumCell];
        cell::dNumSiClxs = new double* [cell::iNumCell];
        cell::dNumSiClxg = new double* [cell::iNumCell];
        cell::dNumNeutral    = new double [cell::iNumCell];
        cell::dNumMask       = new double [cell::iNumCell];
        cell::dNumIon        = new double [cell::iNumCell];
        cell::ElectricForce = new double* [cell::iNumCell];
        cell::iDim           = new int* [cell::iNumCell];  //2-D array int type with size of iNumCell // store the coordinates of each cell
        cell::iID_NBR        = new int* [cell::iNumCell]; // NBR = neighboring cell

        for(i =0; i<cell::iNumCell; i++)        cell::dNumSiClxs[i] = new double [4];
        for(i =0; i<cell::iNumCell; i++)        cell::dNumSiClxg[i] = new double [4];
        for(i =0; i<cell::iNumCell; i++)        cell::ElectricForce[i] = new double [3];
        for(i=0; i<cell::iNumCell; i++)         cell::iDim[i]    = new int [3]; //for each row in iDim, assign a 1-D array with size of 3 to store the (x y z) coordinate
        for(i=0; i<cell::iNumCell; i++)         cell::iID_NBR[i] = new int [27]; // for each row in iID_NBR, assign a 1-D array with size of 27 to store the neighboring coordinates


        for(i=0; i<cell::iNumCell; i++)
        {
                cell::iDim[i][0]       =     i%ix;                     //--x coordinate
                cell::iDim[i][1]       = int(i/ix)%iy;            //--y coordinate
                cell::iDim[i][2]       = int(i/ix/iy);               //--z coordinate
                cell::iStatus[i]       = 0;                                  //--initialization
                cell::dNumCharge[i]    = 0;
                cell::dNumMaterial[i]  = 0;
                cell::dNumNeutral[i]   = 0;
                cell::dNumMask[i]      = 0;
                cell::dNumIon[i]       = 0;
                cell::dNumSiClxs[i][0] = 0;                         //--number of Si(s)
                cell::dNumSiClxs[i][1] = 0;                          //--number of SiCl(s)
                cell::dNumSiClxs[i][2] = 0;                          //--number of SiCl2(s)
                cell::dNumSiClxs[i][3] = 0;                          //--number of SiCl3(s)
                cell::dNumSiClxg[i][0] = 0;                             //--number of Si(g)
                cell::dNumSiClxg[i][1] = 0;                             //--number of SiCl(g)
                cell::dNumSiClxg[i][2] = 0;                             //--number of SiCl2(g)
                cell::dNumSiClxg[i][3] = 0;                             //--number of SiCl3(g)
                cell::ElectricForce[i][0]  =  0.0;
                cell::ElectricForce[i][1]  =  0.0;
                cell::ElectricForce[i][2]  =  0.0;
        }



        cell::setNBR();
}

#pragma acc routine seq
void cell::setNBR()
{
        int i, j, k, itag;
        int ii,jj;
        //                               00     01     02     03     04     05     06     07     08     09     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26
        int ishiftx[27]={  -1,       0,      1,     -1,       0,       1,      -1,      0,       1,      -1,     0,      1,     -1,     0,     1,      -1,     0,      1,     -1,     0,      1,     -1,      0,       1,     -1,      0,       1};
        int ishifty[27]={  -1,     -1,     -1,      0,       0,      0,        1,       1,       1,      -1,    -1,    -1,      0,     0,     0,       1,      1,      1,     -1,    -1,    -1,      0,      0,      0,       1,       1,       1};
        int ishiftz[27]={  -1,     -1,     -1,     -1,     -1,     -1,      -1,     -1,     -1,       0,     0,      0,     0,     0,     0,       0,     0,     0,       1,      1,      1,       1,       1,       1,       1,       1,       1};

        for(jj = 0; jj < cell::iNumCell; jj++)
        {
                for(ii = 0; ii < 27; ii++)
                {
                       i = cell::iDim[jj][0] + ishiftx[ii];
                       j = cell::iDim[jj][1] + ishifty[ii];
                       k = cell::iDim[jj][2] + ishiftz[ii];


                       if(i == cell::iDimSize[0] ) i = 0;  // periodic boundary condition
                       if(j == cell::iDimSize[1] ) j = 0;
                       if(k == cell::iDimSize[2] ) k = 0;
                       if(i == -1                ) i = cell::iDimSize[0] - 1;
                       if(j == -1                ) j = cell::iDimSize[1] - 1;
                       if(k == -1                ) k = cell::iDimSize[2] - 1;
                       itag = i + (j + k * cell::iDimSize[1]) * cell::iDimSize[0];
                       cell::iID_NBR[jj][ii] = itag;
                }
        }

}




#pragma acc routine seq
void cell::setStatus(int iTagCell, int iST, int iNum)
{
        int iVacuumStat = 0;
       int iSubstrateStat = 1, iMaskStat = 2;
        cell::iStatus[iTagCell]       = iST;
        if(iST == iVacuumStat ) // vacuumdvec
        {
                cell::dNumCharge[iTagCell]    = 0;
                cell::dNumMaterial[iTagCell] = 0;
                cell::dNumNeutral[iTagCell]   = 0;
                cell::dNumMask[iTagCell]      = 0;
                cell::dNumIon[iTagCell]       = 0;
        }
        else if(iST == iSubstrateStat ) //material
        {
                cell::dNumCharge[iTagCell]    = 0;
                cell::dNumMaterial[iTagCell]  = iNum;
                cell::dNumSiClxs[iTagCell][0] = cell::dNumMaterial[iTagCell];
                cell::dNumNeutral[iTagCell]   = 0;
                cell::dNumMask[iTagCell]      = 0;
                cell::dNumIon[iTagCell]       = 0;
         }
         else if(iST == iMaskStat ) //mask
         {
                cell::dNumCharge[iTagCell]    = 0;
                cell::dNumMaterial[iTagCell]  = 0;
                cell::dNumNeutral[iTagCell]   = 0;
                cell::dNumMask[iTagCell]      = iNum;
                cell::dNumIon[iTagCell]       = 0;
          }


        //ctor
}



//input: itag, iPos, incident_V; output: norm_surface_N, norm_reflected_V, grazing_angle, incident_angle
#pragma acc routine seq
void surface_normal(int searching_index [][3], int searching_number, int itag, int* iPos, double* incident_V,
                                              int Nx, int Ny, int Nz, int* iStatus, int** iID_NBR, int iVacuumStat, int iSubstrateStat, int iMaskStat,
                                             double* norm_surface_N, double* norm_reflected_V, double* grazing_angle, double* incident_angle ){
        double Surfacesites [searching_number][3];
        int NN_x, NN_y, NN_z, itagNeighbor;
        int NN_x_for_cal_itag, NN_y_for_cal_itag, NN_z_for_cal_itag;
        int X_dir = 0, Y_dir = 1, Z_dir = 2;
        int indexSurfacesites = 0;
        //--determine surface sites
        for (int i = 0; i < searching_number; i++){
                NN_x = iPos[X_dir] + searching_index[i][0];
                NN_y = iPos[Y_dir] + searching_index[i][1];
                NN_z = iPos[Z_dir] + searching_index[i][2];
                NN_x_for_cal_itag = NN_x;
                NN_y_for_cal_itag = NN_y;
                NN_z_for_cal_itag = NN_z;

                //--periodic boundary condition
                while (NN_x_for_cal_itag >= Nx )  NN_x_for_cal_itag -= Nx;
                while (NN_y_for_cal_itag >= Ny )  NN_y_for_cal_itag -= Ny;
                while (NN_z_for_cal_itag >= Nz )  NN_z_for_cal_itag -= Nz;
                while (NN_x_for_cal_itag < 0                )                                NN_x_for_cal_itag += Nx;
                while (NN_y_for_cal_itag < 0                )                                NN_y_for_cal_itag += Ny;
                while (NN_z_for_cal_itag < 0                )                                NN_z_for_cal_itag += Nz;

                //--calculated itag of Neighboring cell
                itagNeighbor = NN_x_for_cal_itag + ( NN_y_for_cal_itag + NN_z_for_cal_itag * Ny)*Nx;


                //--check if the cell property of center cell and neighboring cell same
                if ( iStatus[itagNeighbor] == iSubstrateStat || iStatus[itagNeighbor] == iMaskStat){
                        //--check six faces of neighboring cells are exposed to vacuum
                        if (  iStatus[iID_NBR[itagNeighbor][4]] == iVacuumStat || iStatus[iID_NBR[itagNeighbor][12]] == iVacuumStat ||
                                iStatus[iID_NBR[itagNeighbor][10]] == iVacuumStat || iStatus[iID_NBR[itagNeighbor][14]] == iVacuumStat ||
                                iStatus[iID_NBR[itagNeighbor][16]] == iVacuumStat || iStatus[iID_NBR[itagNeighbor][22]] == iVacuumStat){
                                Surfacesites[indexSurfacesites][0] = double(NN_x);
                                Surfacesites[indexSurfacesites][1] = double(NN_y);
                                Surfacesites[indexSurfacesites][2] = double(NN_z);
                                indexSurfacesites++;
                        }
                }
        }




        //cin.get();
        double xbar, ybar, zbar, xsum, ysum, zsum, numSites, A11, A12, A13, A21, A22, A23, A31, A32, A33;
        xbar = 0; ybar = 0; zbar = 0; xsum = 0; ysum = 0; zsum = 0; numSites = 0;
        A11 = 0; A12 = 0; A13 = 0; A21 = 0; A22 = 0; A23 = 0; A31 = 0; A32 = 0; A33 = 0;

        for(int i =0; i< indexSurfacesites; i++){
                xsum = xsum + Surfacesites[i][X_dir]*Surfacesites[i][3];
                ysum = ysum + Surfacesites[i][Y_dir]*Surfacesites[i][3];
                zsum = zsum + Surfacesites[i][Z_dir]*Surfacesites[i][3];
                numSites++;
        }

        xbar = xsum/numSites;
        ybar = ysum/numSites;
        zbar = zsum/numSites;

        //--calculated surface normal
        for(int i =0; i< indexSurfacesites; i++){
                A11 = A11 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][X_dir] - xbar)*Surfacesites[i][3];
                A22 = A22 + (Surfacesites[i][Y_dir] - ybar)*(Surfacesites[i][Y_dir] - ybar)*Surfacesites[i][3];
                A33 = A33 + (Surfacesites[i][Z_dir] - zbar)*(Surfacesites[i][Z_dir] - zbar)*Surfacesites[i][3];
                A12 = A12 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][Y_dir] - ybar)*Surfacesites[i][3];
                A13 = A13 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][Z_dir] - zbar)*Surfacesites[i][3];
                A23 = A23 + (Surfacesites[i][Y_dir] - ybar)*(Surfacesites[i][Z_dir] - zbar)*Surfacesites[i][3];
        }

        A21 = A12;
        A31 = A13;
        A32 = A23;

        double detA = A11*A22*A33+A12*A23*A31+A13*A21*A32-A31*A22*A13-A23*A32*A11-A12*A21*A33;
        int ranknumber = 3;
        double A [9] = {A11, A12, A13, A21, A22, A23, A31, A32, A33};  //--matrix
        double D [3]; //--eigenvalue
        double V [9];  //--eigenvector
        int it_max = 500;
        int it_num;
        int rot_num;
        jacobi_eigenvalue( ranknumber, A, it_max, V, D,  it_num,  rot_num );
        int minVal_index = 0;
        double minVal = fabs(D[0]);
        for (int i = 0; i < 3; i++){
                if( fabs(D[i]) < minVal){
                        minVal = fabs(D[i]);
                        minVal_index = i;
                }
        }
        //--normalized surface normal vector
        norm_surface_N[X_dir] = V[minVal_index*3+X_dir];
        norm_surface_N[Y_dir] = V[minVal_index*3+Y_dir];
        norm_surface_N[Z_dir] = V[minVal_index*3+Z_dir];

        //--determine the direction of the surface normal
        int poll_positive = 0;
        int poll_negative = 0;
        int poll_x, poll_y, poll_z, itag_poll;

        //cout << "position = " << iPos[X_dir] << " " << iPos[Y_dir] << " " << iPos[Z_dir] << endl;

        for (int i = 1; i <= 3 ; i++){
                poll_x = floor(iPos[X_dir] + 2*i*norm_surface_N[X_dir]);
                poll_y = floor(iPos[Y_dir] + 2*i*norm_surface_N[Y_dir]);
                poll_z = floor(iPos[Z_dir] + 2*i*norm_surface_N[Z_dir]);
                //--periodic boundary condition
                while ( poll_x >= Nx )  poll_x -= Nx;
                while ( poll_y >= Ny )  poll_y -= Ny;
                while ( poll_z >= Nz )  poll_z -= Nz;
                while ( poll_x < 0       )   poll_x += Nx;
                while ( poll_y < 0       )   poll_y += Ny;
                while ( poll_z < 0       )   poll_z += Nz;
                itag_poll = poll_x + (poll_y + poll_z*Ny)*Nx;
                if ( iStatus[itag_poll] == iVacuumStat ){
                        poll_positive++;
                }
        }

        for (int i = 1; i <= 3 ; i++){
                poll_x = floor(iPos[X_dir] - 2*i*norm_surface_N[X_dir]);
                poll_y = floor(iPos[Y_dir] - 2*i*norm_surface_N[Y_dir]);
                poll_z = floor(iPos[Z_dir] - 2*i*norm_surface_N[Z_dir]);
                //--periodic boundary condition
                while ( poll_x >= Nx )  poll_x -= Nx;
                while ( poll_y >= Ny )  poll_y -= Ny;
                while ( poll_z >= Nz )  poll_z -= Nz;
                while ( poll_x < 0       )   poll_x += Nx;
                while ( poll_y < 0       )   poll_y += Ny;
                while ( poll_z < 0       )   poll_z += Nz;
                itag_poll = poll_x + (poll_y + poll_z*Ny)*Nx;
                if ( iStatus[itag_poll] == iVacuumStat ){
                        poll_negative++;
                }
        }


        if (poll_positive > poll_negative){
                 norm_surface_N[X_dir] = norm_surface_N[X_dir];
                 norm_surface_N[Y_dir] = norm_surface_N[Y_dir];
                 norm_surface_N[Z_dir] = norm_surface_N[Z_dir];
        }else if (poll_positive < poll_negative){
                 norm_surface_N[X_dir] = (-1)*norm_surface_N[X_dir];
                 norm_surface_N[Y_dir] = (-1)*norm_surface_N[Y_dir];
                 norm_surface_N[Z_dir] = (-1)*norm_surface_N[Z_dir];
         }else if (poll_positive == poll_negative){
                 /*
                 random_device rd;
                 default_random_engine generator( rd() );  //--random number generator
                 uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
                 */
                 long iRandTag = 0;
                 double rdd = ran3(&iRandTag);

                 if( rdd < 0.5 ){
                         norm_surface_N[X_dir] = norm_surface_N[X_dir];
                         norm_surface_N[Y_dir] = norm_surface_N[Y_dir];
                         norm_surface_N[Z_dir] = norm_surface_N[Z_dir];
                 }else if ( 0.5 <= rdd && rdd < 1.0){
                         norm_surface_N[X_dir] = (-1)*norm_surface_N[X_dir];
                         norm_surface_N[Y_dir] = (-1)*norm_surface_N[Y_dir];
                         norm_surface_N[Z_dir] = (-1)*norm_surface_N[Z_dir];
                  }
         }


        //--calculate the normalized reflected velocity
        double speed = sqrt( pow(incident_V[X_dir], 2.0) + pow(incident_V[Y_dir], 2.0) + pow(incident_V[Z_dir], 2.0)  );
        double norm_incident_V [3];
        double n_dot_vi;
        if ( speed == 0 ){
                norm_incident_V[X_dir] = 0.0;  //--x component of normalized incident velocity vector
                norm_incident_V[Y_dir] = 0.0;  //--y component of normalized incident velocity vector
                norm_incident_V[Z_dir] = 0.0;  //--z component of normalized incident velocity vector
        }else{
                norm_incident_V[X_dir] = incident_V[X_dir]/speed;  //--x component of normalized incident velocity vector
                norm_incident_V[Y_dir] = incident_V[Y_dir]/speed;  //--y component of normalized incident velocity vector
                norm_incident_V[Z_dir] = incident_V[Z_dir]/speed;  //--z component of normalized incident velocity vector
        }

        n_dot_vi = norm_surface_N[X_dir]*norm_incident_V[X_dir]
                              +norm_surface_N[Y_dir]*norm_incident_V[Y_dir]
                              +norm_surface_N[Z_dir]*norm_incident_V[Z_dir];

        norm_reflected_V[X_dir] = norm_incident_V[X_dir] - 2*n_dot_vi*norm_surface_N[X_dir];
        norm_reflected_V[Y_dir] = norm_incident_V[Y_dir] - 2*n_dot_vi*norm_surface_N[Y_dir];
        norm_reflected_V[Z_dir] = norm_incident_V[Z_dir] - 2*n_dot_vi*norm_surface_N[Z_dir];
        double angle_between_surface_normal_and_velocity = acos(n_dot_vi)*180/PI;

        if ( angle_between_surface_normal_and_velocity  >  90 ){
                 *incident_angle = 180 - angle_between_surface_normal_and_velocity;
        }else{
                *incident_angle = angle_between_surface_normal_and_velocity;
        }
        *grazing_angle =  90 - *incident_angle;

}









