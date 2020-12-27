#include "etching.h"
#include "cell.h"
#include "rand.h"

using namespace std;
#pragma acc routine seq
void cell::initial(int ix, int iy, int iz, double Lx, double Ly, double Lz)
{
        cell::iDimSize[X_dir] = ix;
        cell::iDimSize[Y_dir] = iy;
        cell::iDimSize[Z_dir] = iz;
        cell::dDimLength[X_dir] = Lx;
        cell::dDimLength[Y_dir] = Ly;
        cell::dDimLength[Z_dir] = Lz;
        cell::dDimLength_per_cell[X_dir]= Lx/ix;
        cell::dDimLength_per_cell[Y_dir] = Ly/iy;
        cell::dDimLength_per_cell[Z_dir] = Lz/iz;

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



void cell::write_to_vtk(std::string datatype, int cNx, int cNy, int cNz, double delta, int* iStatus, double* dNumMaterial, double* dNumMask,
                                       double** dNumSiClxs, int iNumMaterial, bool PrintSi, bool PrintSiCl, bool PrintSiCl2, bool PrintSiCl3,
                                       std::string directory, std::string append, int file_index ){

        string Filename;
        string file_index_name = to_string(file_index);
        //cout<<file_index_name;
        Filename = directory+file_index_name+append;
        ofstream vtk(Filename);
        int number_of_cells = cNx * cNy * cNz ;
        int pNx = cNx + 1;
        int pNy = cNy + 1;
        int pNz = cNz + 1;
        int number_of_points = pNx * pNy * pNz;
        int iVacuumStat = 0, iSubstrateStat = 1, iMaskStat = 2;
        vtk << "# vtk DataFile Version 3.0" << endl;
        vtk << "Data : feature_profile\n";
        vtk << "ASCII" << endl;

        //--Written in an UNSTRUCTURED_GRID
        if ( datatype == "UNSTRUCTURED_GRID" ) {

                vtk << "DATASET UNSTRUCTURED_GRID" << endl;
                vtk << "POINTS  " << number_of_points << "  float" << endl;
                //--Generate grid points coordinates
                for(int z = 0; z < pNz; z++){
                        for(int y = 0; y < pNy; y++){
                                for(int x = 0; x < pNx; x++){
                                        vtk << x * delta << "\t" << y * delta << "\t" << z * delta <<endl;
                                }
                       }
                }

                vtk << endl;
                vtk << "CELLS" << "\t" << number_of_cells << "\t" << 9*number_of_cells << endl;  //-- Each cell is determined by 9 values

                int index;
                for(int z = 0; z < cNz; z++){
                        for(int y = 0; y < cNy; y++){
                                for(int x = 0; x < cNx; x++){
                                        index = x + ( y + z * pNy) * pNx;

                                        vtk << "8" << "\t" << index                                             << "\t" << index + 1
                                                               << "\t" << (index + pNx)                            << "\t" << (index + pNx + 1)
                                                               << "\t" << (index + pNx * pNy)               << "\t" << (index + pNx * pNy + 1)
                                                               << "\t" << (index + pNx * pNy + pNx) << "\t" << (index + pNx * pNy + pNx + 1) << endl;
                                }
                        }
                 }

                 vtk << endl;
                 vtk << "CELL_TYPES" << "    " << number_of_cells <<endl;

                 for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << "11" << "\t";
                                }
                        }
                 }

                 vtk << endl;
                 vtk << "CELL_DATA" <<  "    " << number_of_cells <<endl;

                 vtk << endl;

                vtk << "SCALARS dNumMaterial float 1" << endl;
                 vtk << "LOOKUP_TABLE default" << endl;
                for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << dNumMaterial[ x + ( y + z * cNy ) * cNx ] << endl;
                                 }
                         }
                 }

                 vtk << "SCALARS dNumMask float 1" << endl;
                 vtk << "LOOKUP_TABLE default" << endl;
                for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << dNumMask[ x + ( y + z * cNy ) * cNx ] << endl;
                                 }
                         }
                 }

                 vtk << "SCALARS iStatus float 1" << endl;
                 vtk << "LOOKUP_TABLE default" << endl;
                 for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << iStatus[ x + ( y + z * cNy ) * cNx ] << endl;
                                 }
                         }
                 }

                if ( PrintSi ){
                        vtk << "SCALARS dNumSis float 1" << endl;
                        vtk << "LOOKUP_TABLE default" << endl;
                        for(int z = 0; z < cNz; z++){
                                for(int y = 0; y < cNy; y++){
                                        for(int x = 0; x < cNx; x++){
                                                int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][0]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                 if (PrintSiCl){
                         vtk << "SCALARS dNumSiCls float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                for(int y = 0; y < cNy; y++){
                                        for(int x = 0; x < cNx; x++){
                                                int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][1]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                if ( PrintSiCl2){
                         vtk << "SCALARS dNumSiCl2s float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                for(int y = 0; y < cNy; y++){
                                        for(int x = 0; x < cNx; x++){
                                                int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][2]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
               }

               if ( PrintSiCl3){
                         vtk << "SCALARS dNumSiCl3s float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                for(int y = 0; y < cNy; y++){
                                        for(int x = 0; x < cNx; x++){
                                                int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][3]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
               }

                 vtk.close();
                 cout<<"The result is written to "<<Filename<<endl;
        }

        //--Written in a STRUCTURED_POINTS
        else if ( datatype == "STRUCTURED_POINTS" ){


                 vtk << "DATASET STRUCTURED_POINTS" << endl;
                 vtk << "DIMENSIONS  " << cNx << "  " << cNy << "  " << cNz<<endl;
                 vtk << "ORIGIN " <<"0  0  0" <<endl;
                 vtk << "SPACING " << delta << "  " << delta << "  " << delta <<endl;
                 vtk << "POINT_DATA" <<  "    " << number_of_cells <<endl;

                 vtk << endl;

                 vtk << "SCALARS dNumMaterial float 1" << endl;
                 vtk << "LOOKUP_TABLE default" << endl;
                  for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << dNumMaterial[ x + ( y + z * cNy ) * cNx ] << endl;
                                 }
                         }
                 }


                 vtk << "SCALARS iStatus float 1" << endl;
                 vtk << "LOOKUP_TABLE default" << endl;
                 for(int z = 0; z < cNz; z++){
                         for(int y = 0; y < cNy; y++){
                                 for(int x = 0; x < cNx; x++){
                                         vtk << iStatus[ x + ( y + z * cNy ) * cNx ] << endl;
                                 }
                         }
                 }

                 if ( PrintSi){
                         vtk << "SCALARS dNumSis float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                 for(int y = 0; y < cNy; y++){
                                         for(int x = 0; x < cNx; x++){
                                                 int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][0]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                 if ( PrintSiCl){
                         vtk << "SCALARS dNumSiCls float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                 for(int y = 0; y < cNy; y++){
                                         for(int x = 0; x < cNx; x++){
                                                 int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][1]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                 if ( PrintSiCl2){
                         vtk << "SCALARS dNumSiCl2s float 1" << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                 for(int y = 0; y < cNy; y++){
                                         for(int x = 0; x < cNx; x++){
                                                 int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                        vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][2]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                 if ( PrintSiCl3){
                         vtk << "SCALARS dNumSiCl3s float 1"  << endl;
                         vtk << "LOOKUP_TABLE default" << endl;
                         for(int z = 0; z < cNz; z++){
                                 for(int y = 0; y < cNy; y++){
                                         for(int x = 0; x < cNx; x++){
                                                 int status = iStatus[ x + ( y + z * cNy ) * cNx ];
                                                 if ( status == iVacuumStat  ){
                                                         vtk << iVacuumStat << endl;
                                                 }else if (  status == iMaskStat){
                                                         vtk << iMaskStat << endl;
                                                 }
                                                 else if ( status == iSubstrateStat ){
                                                         vtk << dNumSiClxs[ x + ( y + z * cNy ) * cNx ][3]/iNumMaterial << endl;
                                                 }
                                         }
                                 }
                         }
                 }

                 vtk.close();
                 cout<<"The result is written to "<<Filename<<endl;
        }

}



vector<double> read_data(string filename){

        vector<double> data_read;
        string buffer_string_line;
        string data;
        ifstream in(filename);
	    while (!in.eof())
        {
                getline(in, data);
                istringstream delim_data_input(  data  );
                if (  data.empty() )   continue;
                while(delim_data_input>>buffer_string_line){
                        data_read.push_back(stod(buffer_string_line));
                }
	    }
	    in.close();
        return data_read;
}




