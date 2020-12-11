#include "etching.h"
#include "cell.h"
#include "jacobi_eigenvalue.h"

using namespace std;
#pragma acc routine seq
void cell::initial(int ix, int iy, int iz, double Lx, double Ly, double Lz)
{
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
void cell::surface_normal(int searching_index[][3], int searching_number, int itag, int* iPos,
                                                        double* incident_V, double* norm_surface_N, double* norm_reflected_V, double* grazing_angle, double* incident_angle ){
        vector<vector<double>> Surfacesites ;
        int NN_x, NN_y, NN_z, itagNeighbor;
        int NN_x_for_cal_itag, NN_y_for_cal_itag, NN_z_for_cal_itag;

        //--determine surface sites
        for (int i = 0; i < searching_number; i++){
                NN_x = iPos[X_dir] + searching_index[i][0];
                NN_y = iPos[Y_dir] + searching_index[i][1];
                NN_z = iPos[Z_dir] + searching_index[i][2];
                NN_x_for_cal_itag = NN_x;
                NN_y_for_cal_itag = NN_y;
                NN_z_for_cal_itag = NN_z;

                //--periodic boundary condition
                while (NN_x_for_cal_itag >= cell::iDimSize[X_dir] )  NN_x_for_cal_itag -= cell::iDimSize[X_dir];
                while (NN_y_for_cal_itag >= cell::iDimSize[Y_dir] )  NN_y_for_cal_itag -= cell::iDimSize[Y_dir];
                while (NN_z_for_cal_itag >= cell::iDimSize[Z_dir] )  NN_z_for_cal_itag -= cell::iDimSize[Z_dir];
                while (NN_x_for_cal_itag < 0                )                                NN_x_for_cal_itag += cell::iDimSize[X_dir];
                while (NN_y_for_cal_itag < 0                )                                NN_y_for_cal_itag += cell::iDimSize[Y_dir];
                while (NN_z_for_cal_itag < 0                )                                NN_z_for_cal_itag += cell::iDimSize[Z_dir];

                //--calculated itag of Neighboring cell
                itagNeighbor = NN_x_for_cal_itag + ( NN_y_for_cal_itag + NN_z_for_cal_itag * cell::iDimSize[1])*cell::iDimSize[0];


                //--check if the cell property of center cell and neighboring cell same
                if (cell::iStatus[itagNeighbor] == iSubstrateStat ){
                        //--check six faces of neighboring cells are exposed to vacuum
                        if (  cell::iStatus[cell::iID_NBR[itagNeighbor][4]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][12]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[itagNeighbor][10]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][14]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[itagNeighbor][16]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][22]] == iVacuumStat){
                                Surfacesites.push_back({double(NN_x), double(NN_y), double(NN_z), cell::dNumMaterial[itagNeighbor]});
                                //cout << "dNumMaterial[itagNeighbor] = " << dNumMaterial[itagNeighbor] << endl;
                                //cin.get();
                        }
                }else if ( cell::iStatus[itagNeighbor] == iMaskStat  ){
                        if (  cell::iStatus[cell::iID_NBR[itagNeighbor][4]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][12]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[itagNeighbor][10]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][14]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[itagNeighbor][16]] == iVacuumStat || cell::iStatus[cell::iID_NBR[itagNeighbor][22]] == iVacuumStat){
                                Surfacesites.push_back({double(NN_x), double(NN_y), double(NN_z), cell::dNumMask[itagNeighbor]});
                                //cout << "dNumMask[itagNeighbor] = " << dNumMask[itagNeighbor] << endl;
                                //cin.get();
                        }
                }
        }




        //cin.get();
        double xbar, ybar, zbar, xsum, ysum, zsum, numSites, A11, A12, A13, A21, A22, A23, A31, A32, A33;
        xbar = 0; ybar = 0; zbar = 0; xsum = 0; ysum = 0; zsum = 0; numSites = 0;
        A11 = 0; A12 = 0; A13 = 0; A21 = 0; A22 = 0; A23 = 0; A31 = 0; A32 = 0; A33 = 0;

        for(int i =0; i< Surfacesites.size(); i++){
                xsum = xsum + Surfacesites[i][X_dir]*Surfacesites[i][3];
                ysum = ysum + Surfacesites[i][Y_dir]*Surfacesites[i][3];
                zsum = zsum + Surfacesites[i][Z_dir]*Surfacesites[i][3];
                numSites = numSites + Surfacesites[i][3];
        }
        if ( numSites == 0)  cout << "numSites = 0" << endl;
        xbar = xsum/numSites;
        ybar = ysum/numSites;
        zbar = zsum/numSites;

        //--calculated surface normal
        for(int i = 0; i < Surfacesites.size(); i++){
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
        /*
        if ( detA != 0){
                cout << "norm_surface_N[X_dir] = " << norm_surface_N[X_dir] << endl;
                cout << "norm_surface_N[Y_dir] = " << norm_surface_N[Y_dir] << endl;
                cout << "norm_surface_N[Z_dir] = " << norm_surface_N[Z_dir] << endl;
                cout << "V[0] V[1] V[2] = " << V[0] << "  " << V[1] << "  " << V[2] << "  " << D[0] << endl;
                cout << "V[3] V[4] V[5] = " << V[3] << "  " << V[4] << "  " << V[5] << "  " << D[1] << endl;
                cout << "V[6] V[7] V[8] = " << V[6] << "  " << V[7] << "  " << V[8] << "  " << D[2] << endl;
                cout << A11 << "    " << A12 << "    " << A13 << endl;
                cout << A21 << "    " << A22 << "    " << A23 << endl;
                cout << A31 << "    " << A32 << "    " << A33 << endl;
                cout << "V'[0] V'[1] V'[2] = " << (A11*V[0]+A12*V[1]+A13*V[2])/V[0] <<  "  " << (A21*V[0]+A22*V[1]+A23*V[2])/V[1] << "  " << (A31*V[0]+A32*V[1]+A33*V[2])/V[2] << endl;
                cout << "V'[3] V'[4] V'[5] = " << (A11*V[3]+A12*V[4]+A13*V[5])/V[3] <<  "  " << (A21*V[3]+A22*V[4]+A23*V[5])/V[4] << "  " << (A31*V[3]+A32*V[4]+A33*V[5])/V[5] << endl;
                cout << "V'[6] V'[7] V'[8] = " << (A11*V[6]+A12*V[7]+A13*V[8])/V[6] <<  "  " << (A21*V[6]+A22*V[7]+A23*V[8])/V[7] << "  " << (A31*V[6]+A32*V[7]+A33*V[8])/V[8] << endl;
                for(int i =0; i< Surfacesites.size(); i++)    cout << Surfacesites[i][X_dir] << " " << Surfacesites[i][Y_dir] << " " << Surfacesites[i][Z_dir] << endl;
                cin.get();
        }
        */
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
                while ( poll_x >= cell::iDimSize[X_dir] )  poll_x -= cell::iDimSize[X_dir];
                while ( poll_y >= cell::iDimSize[Y_dir] )  poll_y -= cell::iDimSize[Y_dir];
                while ( poll_z >= cell::iDimSize[Z_dir] )  poll_z -= cell::iDimSize[Z_dir];
                while ( poll_x < 0                )                                poll_x += cell::iDimSize[X_dir];
                while ( poll_y < 0                )                                poll_y += cell::iDimSize[Y_dir];
                while ( poll_z < 0                )                                poll_z += cell::iDimSize[Z_dir];
                /*
                cout << "poll positive :" << endl;
                cout << poll_x << " " << poll_y << " " << poll_z << endl;
                */
                itag_poll = poll_x + (poll_y + poll_z*cell::iDimSize[1])*cell::iDimSize[0];
                if ( cell::iStatus[itag_poll] == iVacuumStat ){
                        poll_positive++;
                }
        }

        for (int i = 1; i <= 3 ; i++){
                poll_x = floor(iPos[X_dir] - 2*i*norm_surface_N[X_dir]);
                poll_y = floor(iPos[Y_dir] - 2*i*norm_surface_N[Y_dir]);
                poll_z = floor(iPos[Z_dir] - 2*i*norm_surface_N[Z_dir]);
                //--periodic boundary condition
                while ( poll_x >= cell::iDimSize[X_dir] )  poll_x -= cell::iDimSize[X_dir];
                while ( poll_y >= cell::iDimSize[Y_dir] )  poll_y -= cell::iDimSize[Y_dir];
                while ( poll_z >= cell::iDimSize[Z_dir] )  poll_z -= cell::iDimSize[Z_dir];
                while ( poll_x < 0                )                                poll_x += cell::iDimSize[X_dir];
                while ( poll_y < 0                )                                poll_y += cell::iDimSize[Y_dir];
                while ( poll_z < 0                )                                poll_z += cell::iDimSize[Z_dir];
                /*
                cout << "poll negative :" << endl;
                cout << poll_x << " " << poll_y << " " << poll_z << endl;
                */
                itag_poll = poll_x + (poll_y + poll_z*cell::iDimSize[1])*cell::iDimSize[0];
                if ( cell::iStatus[itag_poll] == iVacuumStat ){
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
                 random_device rd;
                 default_random_engine generator( rd() );  //--random number generator
                 uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
                 double rdd = unif(generator);

                 cout << "poll_positive = " << poll_positive << endl;
                 cout << "poll_negative = " << poll_negative << endl;
                 cout << "norm surface N before = " << norm_surface_N[X_dir] << " " << norm_surface_N[Y_dir] << " " << norm_surface_N[Z_dir] << endl;

                 if( rdd < 0.5 ){
                         norm_surface_N[X_dir] = norm_surface_N[X_dir];
                         norm_surface_N[Y_dir] = norm_surface_N[Y_dir];
                         norm_surface_N[Z_dir] = norm_surface_N[Z_dir];
                 }else if ( 0.5 <= rdd && rdd < 1.0){
                         norm_surface_N[X_dir] = (-1)*norm_surface_N[X_dir];
                         norm_surface_N[Y_dir] = (-1)*norm_surface_N[Y_dir];
                         norm_surface_N[Z_dir] = (-1)*norm_surface_N[Z_dir];
                  }

                  //cout << "norm surface N after = " << norm_surface_N[X_dir] << " " << norm_surface_N[Y_dir] << " " << norm_surface_N[Z_dir] << endl;
                  //cin >> poll_negative ;
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

        /*
        if ( cell::iStatus[itag] == iMaskStat){
                for (int i = 0; i < Surfacesites.size(); i++){
                        cout << Surfacesites[i][0] << " " << Surfacesites[i][1] << " " << Surfacesites[i][2] << endl;
                }
                cout << "norm X = " << norm_surface_N[X_dir] <<endl;
                cout << "norm Y = " << norm_surface_N[Y_dir] <<endl;
                cout << "norm Z = " << norm_surface_N[Z_dir] <<endl;
                cout << "incident angle = " << *incident_angle << endl;
                cin.get();
        }
        */


}




void cell::setStatusAll(int iST, int iNum)
{
        int i;
        for(i=0; i<cell::iNumCell; i++)
        {
                cell::setStatus(i,iST,iNum);
        }
}

void cell::setStatusPlane(int inorm, int iConst, int iST, int iNum)
{
        int i[3];
        int iTagCell;
        for(i[0]=0; i[0]<cell::iDimSize[0]; i[0]++)
        {
                for(i[1]=0; i[1]<cell::iDimSize[1]; i[1]++)
                {
                        for(i[2]=0; i[2]<cell::iDimSize[2];i[2]++)
                        {
                                if(i[inorm] == iConst)
                                {
                                        iTagCell = i[0] + (i[1] + i[2] * cell::iDimSize[1]) * cell::iDimSize[0];
                                        cell::setStatus(iTagCell,iST,iNum);
                                }
                        }
                }
        }
}


void cell::setStatusLine(int *idir, int *ipos, int iLength, int iST, int iNum)
{
    int ii,jj;
    int increasement;
    int iTagCell;
        increasement = idir[0] + (idir[1] + idir[2] * cell::iDimSize[1]) * cell::iDimSize[0];
        jj           = ipos[0] + (ipos[1] + ipos[2] * cell::iDimSize[1]) * cell::iDimSize[0];
        for (ii=0; ii<iLength; ii++)
        {
            iTagCell = jj + ii * increasement;
            cell::setStatus(iTagCell,iST, iNum);
        }
}


void cell::setInterFace(int i){
    int j;
        if(cell::iStatus[i] == iSubstrateStat || cell::iStatus[i] == iMaskStat ){
            for (j=0;j<27;j++){
                if(cell::iStatus[cell::iID_NBR[i][j]] == iVacuumStat){
                    cell::iStatus[i] = iInterS_VStat;
                }
            }
        }
}


void cell::setInterFace2(int i){
    int j, iTagNbr;
    int ix[6] = {12, 14, 10 , 16, 4, 22};
        if(cell::iStatus[i] == iVacuumStat){
            for (j=0;j<6;j++){
                iTagNbr = cell::iID_NBR[i][ix[j]];
                if(cell::iStatus[iTagNbr] == iSubstrateStat || cell::iStatus[iTagNbr] == iMaskStat){
                    cell::iStatus[iTagNbr] = iInterS_VStat;
                }
            }
        }
}


void cell::setInterFace4Vacuum(int i){
    int j, iTagNbr;
    int ix[6] = {12, 14, 10 , 16, 4, 22};
    if(cell::iStatus[i] == iVacuumStat){
        for (j=0;j<6;j++){
            iTagNbr = cell::iID_NBR[i][ix[j]];
            if(cell::iStatus[iTagNbr] == iSubstrateStat){
                cell::iStatus[iTagNbr] = iInterS_VStat;
            }else if(cell::iStatus[iTagNbr] == iMaskStat){
                cell::iStatus[iTagNbr] = iInterM_VStat;
            }
        }
    }
}


void cell::TransNeutral4Vacuum(int i){
    int j;
    int ix[6] = {12, 14, 10 , 16, 4, 22};
    int itmpcont = 0;

    for (j=0;j<6;j++){
        if(cell::iStatus[cell::iID_NBR[i][ix[j]]] == iSubstrateStat){
            itmpcont += 1;
        }
    }
    for (j=0;j<6;j++){
        if(cell::iStatus[cell::iID_NBR[i][ix[j]]] == iSubstrateStat){
            cell::dNumNeutral[cell::iID_NBR[i][ix[j]]] = cell::dNumNeutral[i] / itmpcont;
        }
    }
}











