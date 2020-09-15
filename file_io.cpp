#include "cell.h"
#include "etching.h"

using namespace std;


void write_to_vtk(std::string datatype, int cNx, int cNy, int cNz, double delta, int* iStatus, double* dNumMaterial, double* dNumMask,
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



