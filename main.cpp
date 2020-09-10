#include <stdlib.h>
#include <iomanip>
#include <string>
#include "etching.h"
#include "cell.h"
#include "particle.h"
#include "json.hpp"
#include "rand.h"

using namespace std;
using json = nlohmann::json;

double MassChlorine = 35.45*1.66053904E-27;  // unit: kg
double MassArgon = 39.948*1.66053904E-27;  // unit: kg
double MassSilicon = 28.0855*1.66053904E-27;  // unit: kg
double MassElectron = 9.10938356E-31; // unit: kg
int iVacuumStat = 0, iSubstrateStat = 1, iMaskStat = 2, iInterS_VStat = 3, iInterM_VStat = 4;
int P_sputtering = 0;
int C_sputtering = 1;
int X_dir = 0, Y_dir = 1, Z_dir = 2;
const int iClRadicalType = 1, iSigType = 2, iSiClgType = 3, iSiCl2gType = 4, iSiCl3gType = 5, iSiCl4gType = 6;
const int iClIonType = 7, iCl2IonType = 8, iArIonType = 9;
double E0 = 100; //--unit: eV


int main(int argc, char* argv[])
{
        /* Declaration*/
        string inputJsonFile, tmpString, directory, boundary_condition, output_type, append, vtkDataType, geometry, meshfile, meshDATASET ;
        string ion_energy_data_file, ion_angle_data_file, cumulativeflux_ArIon_file, cumulativeflux_ClIon_file, cumulativeflux_Cl2Ion_file;
        string output_file_name, IEADF_data_directory;
        int TotalParticle, timestep_number, output_file_number, Nx, Ny, Nz, searching_radius, searching_number ;
        int iSubstrateThickZ, iMaskThickZ, iTrenchWidthX,  iMaskWidthX, iNumMaterial, iNumMask, particleNumber, file_index, frequency ;
        bool ION_THETA_GAUSSIAN, PrintSi, PrintSiCl, PrintSiCl2, PrintSiCl3 ;
        double dx, dy, dz, Lx, Ly, Lz, dSubstrateThickZ, dMaskThickZ, dTrenchWidthX, dMaskWidthX, Temperature ;
        vector<double> cumulativeflux_ArIon, cumulativeflux_ClIon, cumulativeflux_Cl2Ion, ion_energy, ion_angle;

         /* Initialization*/
        int ThermalParticleTypes = 6; //--Total number of thermal species, here we include Cl radical, Si(g), SiCl(g), SiCl2(g), and SiCl3(g)
        //                                                                                                                          Cl                     Si                    SiCl               SiCl2              SiCl3                    Cl+                Cl2+              Ar+
        //vector<double> ParticleProb_for_inciden_particle = {0.92085,     0.00000,    0.00000,      0.00000,     0.00000,         0.04186,     0.03474,     0.00255};
        vector<double> ParticleProb_for_incident_particle = {0.95238,         0.0,                 0.0,                 0.0,                 0.0,                0.02518,     0.020900,  0.001534};
        double speed_cutoff_for_thermal_paricle [ThermalParticleTypes]; //--Array to store pre-calculated cutoff speed
        vector<double> mass_thermal_particle = { 0.0, MassChlorine, MassSilicon, MassSilicon+MassChlorine, MassSilicon+2*MassChlorine,
                                                                                                     MassSilicon+3*MassChlorine };
        vector<double> phys_sputter_prob =  {0.55,  0.555,   0.60,   0.555,   0.85,    1.3,   1.355,    1.0,   0.75,   0.0};
        vector<double> chem_sputter_prob = {1.00,  1.000,   1.00,   1.000,    1.00,   0.9,  0.800,   0.6,   0.30,   0.0};
        vector<double> p0_ClRadicalReaction = {0.99, 0.40, 0.30, 0.02, 0.0001, 0.08};
        vector<double> p0_redeposition =         {0.02, 0.02, 0.02};
        vector<double> p0_ClIonReaction =      {0.05, 0.10, 0.20, 0.50, 0.50};
        vector<double> p0_Cl2IonReaction =   {0.02, 0.20, 0.25, 0.25, 0.25, 0.25};
        vector<double> p0_ArIonReaction =     {0.05, 0.20, 0.50, 0.50};
        vector<double> Eth_ClIonReaction =    {25, 35, 10, 10, 10};
        vector<double> Eth_Cl2IonReaction = {25, 10, 10, 10, 10, 10};
        vector<double> Eth_ArIonReaction =   {25, 10, 10, 10};
        vector<int> type_ClIonReaction = {P_sputtering, P_sputtering, C_sputtering, C_sputtering, C_sputtering};
        vector<int> type_Cl2IonReaction = {P_sputtering, C_sputtering, C_sputtering, C_sputtering, C_sputtering, C_sputtering};
        vector<int> type_ArIonReaction = {P_sputtering, C_sputtering, C_sputtering, C_sputtering};
        double E0_ClIonReaction = 100;
        double E0_Cl2IonReaction = 100;
        double E0_ArIonReaction = 100;
        int reaction_number[26] = {0};  //--for reaction index from 1 - 26
        int n_for_cosine_law = 1;
        double particle_radius_scale = 0.2;

        /*Input Json file*/
        inputJsonFile = "none" ;
        if ( argc > 0 ){
                tmpString = argv[1] ;
                if ( tmpString == "-input_json" ){
	                    if ( argc > 2 ){
	                            cout << "Input file name found." << endl ;
	                            tmpString = argv[2] ;
		                        if ( tmpString.find( ".json", 0 ) != string::npos ){
	                                    cout << "OK. File name with .json" << endl ;
		                                inputJsonFile = argv[2] ;
		                        }else{
                                        cout << "Error. File name without .json" << endl ;
	                            }
	                    }else{
		                        cout << "Input file name not found." << endl ;
	                    }
                }else{
                        cout << "-input_json tag not found." << endl ;
                }
        }
        cout << "Input Json File Name : " << inputJsonFile << endl ;

        if ( inputJsonFile == "none" ){
                cout << "Please check the input file." << endl;
        }else{   //--Read Json file
                json configuration_json ; //--Declare a json object
                ifstream jsonFile ;
                //jsonFile.open( inputJsonFile.c_str(), ifstream::in ) ; //--read json file named inputJsonFile
                jsonFile.open( inputJsonFile, ifstream::in ) ; //--read json file named inputJsonFile
                jsonFile >> configuration_json ;

                //--mesh setting
                if ( configuration_json.find( "mesh" ) != configuration_json.end() ){
                        json config = configuration_json[ "mesh" ] ;
                        if ( config.find( "geometry" ) != config.end() ){
                                geometry = config[ "geometry" ] ;
	                    }else{
	                            cout << "No geometry tag in input file! " << endl ;
	                            exit( 1 ) ;
	                    }
	                    if ( config.find( "meshfile" ) != config.end() ){
	                            meshfile = config[ "meshfile" ] ;
	                    }else{
	                            cout << "No meshfile tag in input file! " << endl ;
	                            //meshfile = "none" ;
		                        exit( 1 ) ;
	                    }
	                    dx                                 = config["cell_size"] ; //--unit : m, x length per cell size
	                    Nx                                = config["domain_size_x"]; //--number of cells in x direction
                        Ny                                = config["domain_size_y"]; //--number of cells in y direction
                        Nz                                = config["domain_size_z"]; //--number of cells in z direction
                        iSubstrateThickZ = config["substrate_thickness"];
                        iMaskThickZ           = config["mask_thickness"];
                        iTrenchWidthX      = config["trench_width"];
                        iMaskWidthX         = config["mask_width"];
                        iNumMaterial        = config["atom_number_in_material"];
                        iNumMask              = config["atom_number_in_mask"];
                        searching_radius = config["searching_radius"] ;
                }else{
	                    cout << "No mesh tag in input file! " << endl ;
	                    exit( 1 ) ;
                }

                //--Output setting
                if ( configuration_json.find( "output_setting" ) != configuration_json.end() ){
                        json config = configuration_json[ "output_setting" ] ;
                        if ( config[ "type_name" ] == "TECPLOT" ){
                                output_type =	"TECPLOT" ;
                        }else if ( config[ "type_name" ] == "VTK" ){
	                            output_type	= "VTK" ;
	                            append = ".vtk";
	                    }else{
                                output_type	= "VTK" ;
                                append = ".vtk";
	                    }
	                    vtkDataType = config["vtkDataType"] ;
	                    directory = config["output_file_directory"] ;
	                    output_file_name = config["output_file_name"];
                        output_file_number = config["output_file_number"] ;
                }else{
	                    output_type = "VTK" ;
                }

                //--simulation conditions setting
                if ( configuration_json.find( "simulation_conditions" ) != configuration_json.end() ){
                        json config = configuration_json[ "simulation_conditions" ] ;
                        //cout << config << endl;

                        if ( config.find( "timestep_number" ) != config.end() ){
	                            timestep_number = config[ "timestep_number" ] ;
                        }
	                    TotalParticle = config["TotalParticle"];

                        if ( config.find( "boundary_condition" ) != config.end() ){
                                boundary_condition = config["boundary_condition"];
	                    }else{
	                            boundary_condition = "periodic" ;
	                    }
	                    Temperature = config["Temperature"];
                        IEADF_data_directory = config["IEADF_data_directory"];
                        ION_THETA_GAUSSIAN = config["ION_THETA_GAUSSIAN"];
                        PrintSi = config["PrintSi"];
                        PrintSiCl = config["PrintSiCl"];
                        PrintSiCl2 = config["PrintSiCl2"];
                        PrintSiCl3 = config["PrintSiCl3"];
                }
                cout << endl;
                cout << "JSON FILE OUTPUT : " << endl;
                cout << "Mesh:" << endl ;
                cout << "  geometry	               =    " << geometry << endl ;
                cout << "  meshfile	               =    " << meshfile << endl ;
                cout << "  cell_size                    =    " << dx << endl ;
                cout << "  domain size x                =    " << Nx << endl ;
                cout << "  domain size y                =    " << Ny << endl ;
                cout << "  domain size z                =    " << Nz << endl ;
                cout << "  substrate_thickness          =    " << iSubstrateThickZ << endl;
                cout << "  mask_thickness               =    " << iMaskThickZ << endl;
                cout << "  trench_width                 =    " << iTrenchWidthX << endl;
                cout << "  mask_width                   =    " << iMaskWidthX << endl;
                cout << "  atom_number_in_material      =    " << iNumMaterial << endl;
                cout << "  atom_number_in_mask          =    " << iNumMask << endl;
                cout << "  searching_radius             =    " << searching_radius << endl;
                cout << endl;
                cout << "Output:" << endl ;
                cout << "  form                         =    " << output_type << endl ;
                cout << "  vtkDataType                  =    " << vtkDataType << endl ;
                cout <<"  output_file_directory        =    "<< directory << endl ;
                cout <<"  output_file_name             =    "<< output_file_name << endl ;
                cout << "  output_file_number           =    "<< output_file_number << endl;
                cout << endl;
                cout << "Simulation conditions:" << endl ;
                cout << "  TotalParticle                =    " << TotalParticle << endl ;
                cout << "  timestep_number              =    " << timestep_number << endl ;
                cout << "  Boundary condition           =    " << boundary_condition << endl;
                cout << "  Temperature                  =    "  << Temperature << endl;
                cout << "  IEADF_data_directory         =    " << IEADF_data_directory << endl;
                cout << "  ION_THETA_GAUSSIAN           =    " << ION_THETA_GAUSSIAN << endl;
                cout << "  PrintSi                      =    " << PrintSi << endl ;
                cout << "  PrintSiCl                    =    " << PrintSiCl << endl ;
                cout << "  PrintSiCl2                   =    " << PrintSiCl2 << endl ;
                cout << "  PrintSiCl3                   =    " << PrintSiCl3 << endl ;
                cout << endl;
        }


        /*Read IEADF data*/
        ion_energy_data_file =            IEADF_data_directory + "energy.dat";
        ion_angle_data_file =               IEADF_data_directory + "angle.dat";
        cumulativeflux_ClIon_file =   IEADF_data_directory + "cumulativefluxClIon_data.dat";
        cumulativeflux_Cl2Ion_file = IEADF_data_directory + "cumulativefluxCl2Ion_data.dat";
        cumulativeflux_ArIon_file =   IEADF_data_directory + "cumulativefluxArIon_data.dat";
        ion_energy                         = read_data(   ion_energy_data_file    );
        ion_angle                            = read_data(   ion_angle_data_file    );
        cumulativeflux_ClIon    = read_data(   cumulativeflux_ClIon_file   );
        cumulativeflux_Cl2Ion = read_data(   cumulativeflux_Cl2Ion_file  );
        cumulativeflux_ArIon   = read_data(   cumulativeflux_ArIon_file   );

        /*Initialization for file output*/
        file_index = 0;
        frequency = TotalParticle/output_file_number;
        particleNumber = 0;

        /* Pre-calculation of speed cutoff and time interval for each particle type */
        speed_cutoff_for_thermal_paricle[iClRadicalType] = calc_v_cut(  mass_thermal_particle.at(iClRadicalType),   Temperature  );
        speed_cutoff_for_thermal_paricle[iSigType]              = calc_v_cut(  mass_thermal_particle.at(iSigType),                Temperature  );
        speed_cutoff_for_thermal_paricle[iSiClgType]          = calc_v_cut(  mass_thermal_particle.at(iSiClgType),           Temperature  );
        speed_cutoff_for_thermal_paricle[iSiCl2gType]       = calc_v_cut(  mass_thermal_particle.at(iSiCl2gType),         Temperature  );
        speed_cutoff_for_thermal_paricle[iSiCl3gType]       = calc_v_cut(  mass_thermal_particle.at(iSiCl3gType),         Temperature  );

        /*Pre-calculation for surface sites*/
        int searching_range = 2*searching_radius+1;
        searching_number = pow(  searching_range, 3);
        //cout << searching_number << endl;
        //cin.get();
        int searching_index [searching_number][3];
        for (int k = 0; k < searching_range; k++){
                for (int j = 0; j < searching_range; j++){
                        for(int i = 0; i < searching_range; i++){
                                searching_index[i+(j+k*searching_range )*searching_range ][0] = i - searching_radius;
                                searching_index[i+(j+k*searching_range )*searching_range ][1] = j - searching_radius;
                                searching_index[i+(j+k*searching_range )*searching_range ][2] = k - searching_radius;
                        }
                }
        }

        /* Initialization for cell */
        cell C1;
        if ( meshfile == "none" ){
                dy = dx; //--y length per cell size
                dz = dx; //--z length per cell size
                Lx = Nx*dx;
                Ly = Ny*dy;
                Lz = Nz*dz;
                dSubstrateThickZ = iSubstrateThickZ*dz;
                dMaskThickZ = iMaskThickZ*dz;
                dTrenchWidthX = iTrenchWidthX*dx;
                dMaskWidthX = iMaskWidthX*dx;
                C1.initial(Nx, Ny, Nz, Lx, Ly, Lz);


                for(int iz = 0; iz < Nz; iz++){
                        for(int iy = 0; iy < Ny; iy++){
                                for(int ix = 0; ix < Nx; ix++){
                                        int itag =  ix + ( iy + iz*Ny )*Nx;
                                        if(  iz < iSubstrateThickZ ){
                                                //--setting for substrate
                                                C1.setStatus(itag, iSubstrateStat, iNumMaterial);     //--Si has 8 atoms per unit cell with the dimension of 0.54 nm
                                        }else if (  iz >= iSubstrateThickZ && iz < (iSubstrateThickZ + iMaskThickZ)  ) {
                                                if(   ix > iMaskWidthX/2 && ix < (iMaskWidthX/2+iTrenchWidthX)  )   {
                                                        //--setting for trench
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else{
                                                        //--setting for mask
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }
                                        }
                                }
                        }
                }

        }else{
                string buffer_string_line;
                string data;
                ifstream in(meshfile);
                int POINT_DATA, status;
                vector<int> array_NumMaterial;

                if(in.fail()){
                        cout << meshfile << " doest not exist." << endl;
                        exit(1);
                }

	            while (!in.eof())
                {
                        getline(in, data);
                        if (data == "ASCII"){
                                getline(in, data);
                                istringstream delim_data_input(  data  );
                                delim_data_input >> buffer_string_line;
                                if( buffer_string_line == "DATASET")    delim_data_input >> buffer_string_line;
                                meshDATASET =  buffer_string_line;
                                cout << "mesh DATASET = " << meshDATASET<< endl;
                                delim_data_input.clear();
                                getline(in, data);
                                delim_data_input.str(  data  );
                                delim_data_input >> buffer_string_line;
                                if( buffer_string_line == "DIMENSIONS"){
                                        delim_data_input>>buffer_string_line;
                                        Nx = stoi(buffer_string_line) ;
                                        delim_data_input>>buffer_string_line;
                                        Ny = stoi(buffer_string_line)  ;
                                        delim_data_input>>buffer_string_line;
                                        Nz = stoi(buffer_string_line)  ;
                                        delim_data_input.clear();
                                }
                                getline(in, data);
                                getline(in, data);
                                delim_data_input.str(  data  );
                                delim_data_input >> buffer_string_line;
                                if( buffer_string_line == "SPACING"){
                                        delim_data_input>>buffer_string_line;
                                        dx = stod(buffer_string_line) ;
                                        delim_data_input>>buffer_string_line;
                                        dy = stod(buffer_string_line)  ;
                                        delim_data_input>>buffer_string_line;
                                        dz = stod(buffer_string_line)  ;
                                        delim_data_input.clear();
                                }
                                getline(in, data);
                                delim_data_input.str(  data  );
                                delim_data_input >> buffer_string_line;
                                if (buffer_string_line == "POINT_DATA" ){
                                        delim_data_input>>buffer_string_line;
                                        POINT_DATA = stoi(buffer_string_line) ;
                                        delim_data_input.clear();
                                }
                                Lx = Nx*dx;
                                Ly = Ny*dy;
                                Lz = Nz*dz;
                                C1.initial(Nx, Ny, Nz, Lx, Ly, Lz);
                        }

                        if (data == "SCALARS dNumMaterial float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                array_NumMaterial.push_back(stoi(data));
                                                //cout << "data = " << stoi(data) <<endl;
                                                //cin.get();
                                        }
                                }
                         }

                         if (data == "SCALARS iStatus float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                status = stoi(data);
                                                C1.setStatus(i, status, array_NumMaterial[i]);
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSis float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][0] = stod(data);
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSiCls float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][1] = stod(data);
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSiCl2s float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][2] = stod(data);
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSiCl3s float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][3] = stod(data);
                                        }
                                }
                         }
	            }
	            in.close();
        }
        cout << endl;
        cout << "GENERAL INFORMATION : " << endl;
        cout << "meshfile is read from " << meshfile << endl ;
        cout << "Nx = "<<Nx<<", Lx = " << Lx <<endl;
        cout << "Ny = "<<Ny<<", Ly = " << Ly <<endl;
        cout << "Nz = "<<Nz<<", Lz = " << Lz <<endl;
        cout << "speed cutoff for Cl radical = " << speed_cutoff_for_thermal_paricle[iClRadicalType] <<" m/s"<<endl;
        cout << "speed cutoff for Si(g) = "           << speed_cutoff_for_thermal_paricle[iSigType]              <<" m/s"<<endl;
        cout << "speed cutoff for SiCl(g) = "      << speed_cutoff_for_thermal_paricle[iSiClgType]          <<" m/s"<<endl;
        cout << "speed cutoff for SiCl2(g) = "   << speed_cutoff_for_thermal_paricle[iSiCl2gType]        <<" m/s"<<endl;
        cout << "speed cutoff for SiCl3(g) = "   << speed_cutoff_for_thermal_paricle[iSiCl3gType]        <<" m/s"<<endl;
        cout << endl;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumSiClxs, iNumMaterial, PrintSi, PrintSiCl, PrintSiCl2, PrintSiCl3,
                                                               directory+output_file_name, append, file_index );


        //--random number generation
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution


        /*particle size*/
        double dx_particle = particle_radius_scale*dx;
        double dy_particle = particle_radius_scale*dy;
        double dz_particle = particle_radius_scale*dz;
        double six_point_particle [6][3] = {  {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1} };


        #pragma omp parallel for
        for (int indexParticle = 0;   indexParticle < TotalParticle;   indexParticle++){

                particle P1;                                                                                             //--Generate a particle
                P1.setInitialType(ParticleProb_for_incident_particle);  //--choose a particular type of particle
                double rdd1, rdd2, rdd3;                                                                  //--random number for selecting ion energy and angle
                int ReactionExecution = 0;                                                            //--determine if a reaction happen
                int EmitParticle = 0;                                                                          //--index for emitted particle such as SiClx(g)
                int ReflectedParticle = 0;                                                                //--index for original reflected particle such as Ar*, Cl*
                int reaction_index = 0;                                                                    //--index for a particular reaction formula
                double norm_surface_N [3];                                                       //--normalized surface normal vector (x, y, z)
                double norm_reflected_V [3];                                                     //--normalized reflected velocity vector (x, y, z)
                double reflected_velocity [3];                                                      //--reflected velocity (x, y, z)
                double grazing_angle = 0;                                                            //--angle between surface and velocity
                double incident_angle = 0;                                                           //--angle between normal and velocity
                int old_iPos [3];                                                                                   //--to record the iPos (x, y, z) before propogation
                double old_dPos [3];                                                                        //--to record the dPos (x, y, z) before propagation
                int itag;                                                                                                     //--itag of particle center
                int old_itag;
                int itag_six_point [6] ;                                                                     //--itag of particle's six neighboring points
                double dPos_six_point [6][3] ;                                                   //--dPos (x, y, z) of particle's six neighboring points
                int iPos_six_point [6][3] ;                                                              //--iPos (x, y, z) of particle's siz neighboring points
                int count_point_on_solid;                                                             //--to count how many point of a seven-point molecule is on solid cell
                double norm_EmitParticle_V [3];                                             //--normalized velocity for emitted particle
                double dPos_EmitParticle [3];                                                    //--dPos for emitted particle
                int iPos_EmitParticle [3];                                                               //--iPos for emitted particle
                double alpha;                                                                                        //--azimusal angle of surface normal vector
                double beta;                                                                                           //--polar angle of surface normal vector;

                //--Count total number of particle and write to a file in a vtk format
                #pragma omp critical
                {
                        particleNumber++;
                        if(  particleNumber%frequency == 0  ){
                                file_index++;
                                if (file_index == output_file_number){
                                        PrintSi = true;
                                        PrintSiCl = true;
                                        PrintSiCl2 = true;
                                        PrintSiCl3 = true;
                                }
                                write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumSiClxs, iNumMaterial,
                                                                PrintSi, PrintSiCl, PrintSiCl2, PrintSiCl3, directory+output_file_name, append, file_index );
                        }
                }

                if (P1.ParticleType == 0)    continue;

                switch(P1.ParticleType){
                        case iClRadicalType:
                                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                                P1.mass = MassChlorine;
                                P1.speed = P1.setInitialSpeed( Temperature,   MassChlorine,   speed_cutoff_for_thermal_paricle[iClRadicalType]  );
                                P1.energy = 0.5 * P1.mass * P1.speed * P1.speed; //--unit: Joule
                                P1.theta = P1.setInitialTheta(P1.ParticleType); //--unit: radian
                                P1.phi = unif(generator)*2*PI; //--unit: radian
                                break;

                        case iClIonType:
                                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                                P1.mass = MassChlorine;
                                rdd1 = unif(generator);
                                rdd2 = unif(generator);
                                rdd3 = unif(generator);
                                if ( ION_THETA_GAUSSIAN == true){
                                        P1.theta = P1.setInitialTheta(P1.ParticleType); //--unit: radian
                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                }
                                for(int i = 0;    i < cumulativeflux_ClIon.size()-1;     i++){
                                        if  (  rdd1 >= cumulativeflux_ClIon[i] && rdd1 < cumulativeflux_ClIon[i+1]  ){
                                                if ( ION_THETA_GAUSSIAN == false){
                                                        P1.theta =(  ion_angle[i] + rdd2*(ion_angle[i+1] - ion_angle[i])  )*PI/180+PI; //--unit: radian
                                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                                }
                                                P1.energy = (  ion_energy[i] + rdd3*(ion_energy[i+1] - ion_energy[i]) )/Joule_to_eV;//--unit: Joule
                                                P1.speed = sqrt(2*P1.energy/P1.mass);
                                                break;
                                        }
                                }
                                break;

                        case iCl2IonType:
                                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                                P1.mass = MassChlorine*2;
                                rdd1 = unif(generator);
                                rdd2 = unif(generator);
                                rdd3 = unif(generator);
                                if ( ION_THETA_GAUSSIAN == true){
                                        P1.theta = P1.setInitialTheta(P1.ParticleType); //--unit: radian
                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                }
                                for(int i = 0;    i < cumulativeflux_Cl2Ion.size()-1;     i++){
                                        if(  rdd1 >= cumulativeflux_Cl2Ion[i] && rdd1 <  cumulativeflux_Cl2Ion[i+1]  ){
                                                if ( ION_THETA_GAUSSIAN == false){
                                                        P1.theta =(  ion_angle[i] + rdd2*(ion_angle[i+1] - ion_angle[i])  )*PI/180+PI; //--unit: radian
                                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                                }
                                                P1.energy = (  ion_energy[i] + rdd3*(ion_energy[i+1] - ion_energy[i]) )/Joule_to_eV;//--unit: Joule
                                                P1.speed = sqrt(2*P1.energy/P1.mass);
                                                break;
                                        }
                                }
                                break;

                        case  iArIonType:
                                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                                P1.mass = MassArgon;
                                rdd1 = unif(generator);
                                rdd2 = unif(generator);
                                rdd3 = unif(generator);
                                if ( ION_THETA_GAUSSIAN == true){
                                        P1.theta = P1.setInitialTheta(P1.ParticleType); //--unit: radian
                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                }
                                for(int i = 0;  i <  cumulativeflux_ArIon.size()-1;  i++){
                                        if  (  rdd1 >= cumulativeflux_ArIon[i] && rdd1 < cumulativeflux_ArIon[i+1]  ){
                                                if ( ION_THETA_GAUSSIAN == false){
                                                        P1.theta =(  ion_angle[i] + rdd2*(ion_angle[i+1] - ion_angle[i])  )*PI/180+PI; //--unit: radian
                                                        P1.phi = unif(generator)*2*PI; //--unit: radian
                                                }
                                                P1.energy = (  ion_energy[i] + rdd3*(ion_energy[i+1] - ion_energy[i]) )/Joule_to_eV;//--unit: Joule
                                                P1.speed = sqrt(2*P1.energy/P1.mass);
                                                break;
                                        }
                                }
                                break;
                }

                P1.time_interval = dx/P1.speed;
                P1.Vel[X_dir] = P1.speed*sin(P1.theta)*cos(P1.phi);
                P1.Vel[Y_dir] = P1.speed*sin(P1.theta)*sin(P1.phi);
                P1.Vel[Z_dir] = P1.speed*cos(P1.theta);


                for(int indexTimeStep = 0; indexTimeStep < timestep_number; indexTimeStep++){


                        //--chech if particle has been deactivated
                        if(P1.ParticleType == 0){
                                if (EmitParticle == 0){
                                        break;
                                }else{
                                        P1.ParticleType = EmitParticle;
                                        if (EmitParticle == iSiClgType){
                                                P1.speed = P1.setInitialSpeed( Temperature,   MassSilicon+MassChlorine,   speed_cutoff_for_thermal_paricle[iSiClgType]  );
                                        }else if ( EmitParticle == iSiCl2gType){
                                                P1.speed = P1.setInitialSpeed( Temperature,   MassSilicon+2*MassChlorine,   speed_cutoff_for_thermal_paricle[iSiCl2gType]  );
                                        }else if ( EmitParticle == iSiCl3gType){
                                                P1.speed = P1.setInitialSpeed( Temperature,   MassSilicon+3*MassChlorine,   speed_cutoff_for_thermal_paricle[iSiCl3gType]  );
                                        }
                                        P1.dPos[X_dir] = dPos_EmitParticle[X_dir];
                                        P1.dPos[Y_dir] = dPos_EmitParticle[Y_dir];
                                        P1.dPos[Z_dir] = dPos_EmitParticle[Z_dir];
                                        P1.iPos[X_dir] = iPos_EmitParticle[X_dir];
                                        P1.iPos[Y_dir] = iPos_EmitParticle[Y_dir];
                                        P1.iPos[Z_dir] = iPos_EmitParticle[Z_dir];
                                        P1.theta = P1.setReemitTheta(n_for_cosine_law); //--theta is with respect to surface normal
                                        P1.phi = unif(generator)*2*PI;
                                        if( norm_EmitParticle_V[X_dir] == 0 && norm_EmitParticle_V[Y_dir] == 0 ){
                                                alpha = 0;
                                        }else{
                                                alpha = acos(norm_EmitParticle_V[X_dir] / sqrt( pow(norm_EmitParticle_V[X_dir],2.0)+pow(norm_EmitParticle_V[Y_dir],2.0)  ) );
                                        }
                                        beta = acos(norm_EmitParticle_V[Z_dir]);
                                        P1.Vel[X_dir] = P1.speed*( cos(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+cos(beta)*sin(alpha)*cos(P1.theta)-sin(beta)*sin(P1.theta)*sin(P1.phi));
                                        P1.Vel[Y_dir] = P1.speed*( sin(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)+cos(beta)*sin(P1.theta)*sin(P1.phi));
                                        P1.Vel[Z_dir] = P1.speed*( -sin(alpha)*sin(P1.theta)*cos(P1.phi)+cos(alpha)*cos(P1.theta));
                                        EmitParticle = 0;
                                }
                        }


                        //--Give old_iPos and old_dPos values
                        old_iPos[X_dir] = P1.iPos[X_dir];
                        old_iPos[Y_dir] = P1.iPos[Y_dir];
                        old_iPos[Z_dir] = P1.iPos[Z_dir];
                        old_dPos[X_dir] = P1.dPos[X_dir] ;
                        old_dPos[Y_dir] = P1.dPos[Y_dir] ;
                        old_dPos[Z_dir] = P1.dPos[Z_dir] ;
                        old_itag = old_iPos[X_dir] + ( old_iPos[Y_dir] + old_iPos[Z_dir]*C1.iDimSize[Y_dir] )*C1.iDimSize[X_dir];



                        //--propagation
                        P1.Vel[X_dir]       = P1.Vel[X_dir] - C1.grad_potential[old_itag][X_dir] / P1.mass * P1.time_interval;
                        P1.Vel[Y_dir]       = P1.Vel[Y_dir] - C1.grad_potential[old_itag][Y_dir] / P1.mass * P1.time_interval;
                        P1.Vel[Z_dir]       = P1.Vel[Z_dir] - C1.grad_potential[old_itag][Z_dir] / P1.mass * P1.time_interval;
                        P1.dPos[X_dir]   = P1.dPos[X_dir] + P1.Vel[X_dir] * P1.time_interval - 0.5 * C1.grad_potential[old_itag][X_dir] / P1.mass * pow(P1.time_interval, 2);
                        P1.dPos[Y_dir]   = P1.dPos[Y_dir] + P1.Vel[Y_dir] * P1.time_interval - 0.5 * C1.grad_potential[old_itag][Y_dir] / P1.mass * pow(P1.time_interval, 2);
                        P1.dPos[Z_dir]   = P1.dPos[Z_dir] + P1.Vel[Z_dir] * P1.time_interval - 0.5 * C1.grad_potential[old_itag][Z_dir] / P1.mass * pow(P1.time_interval, 2);
                        P1.iPos[X_dir]     = int(   floor(P1.dPos[X_dir]/C1.dDimLength_per_cell[X_dir] )  );
                        P1.iPos[Y_dir]     = int(   floor(P1.dPos[Y_dir] /C1.dDimLength_per_cell[Y_dir] ) );
                        P1.iPos[Z_dir]     = int(   floor(P1.dPos[Z_dir]/C1.dDimLength_per_cell[Z_dir] ) );



                        //--boundary condition
                        if (boundary_condition == "reflective"){
                                if(  P1.iPos[X_dir] <  0 ){
                                        P1.dPos[X_dir] = P1.dPos[X_dir]*(-1);
                                        P1.iPos[X_dir] = floor(P1.dPos[X_dir]/C1.dDimLength_per_cell[X_dir]);
                                        P1.Vel[X_dir] = P1.Vel[X_dir]*(-1);
                                }else if( P1.iPos[X_dir] >= C1.iDimSize[X_dir] ){
                                        P1.dPos[X_dir] = 2*C1.dDimLength[X_dir] - P1.dPos[X_dir];
                                        P1.iPos[X_dir] = floor(P1.dPos[X_dir]/C1.dDimLength_per_cell[X_dir]);
                                        P1.Vel[X_dir] = P1.Vel[X_dir]*(-1);
                                }
                                if( P1.iPos[Y_dir] < 0 ){
                                        P1.dPos[Y_dir] = P1.dPos[Y_dir]*(-1);
                                        P1.iPos[Y_dir] = floor(P1.dPos[Y_dir] /C1.dDimLength_per_cell[Y_dir]);
                                        P1.Vel[Y_dir] = P1.Vel[Y_dir]*(-1);
                                }else if( P1.iPos[Y_dir] >= C1.iDimSize[Y_dir] ){
                                        P1.dPos[Y_dir] = 2*C1.dDimLength[Y_dir] - P1.dPos[Y_dir];
                                        P1.iPos[Y_dir] = floor(P1.dPos[Y_dir] /C1.dDimLength_per_cell[Y_dir]);
                                        P1.Vel[Y_dir] = P1.Vel[Y_dir]*(-1);
                                }
                        }else if (boundary_condition == "periodic"){
                                if(P1.iPos[X_dir] < 0){
                                        P1.dPos[X_dir] = P1.dPos[X_dir] + C1.dDimLength[X_dir];
                                        P1.iPos[X_dir] = P1.iPos[X_dir] + C1.iDimSize[X_dir];
                                }else if(P1.iPos[X_dir] >= C1.iDimSize[X_dir]){
                                        P1.dPos[X_dir] = P1.dPos[X_dir] - C1.dDimLength[X_dir];
                                        P1.iPos[X_dir] = P1.iPos[X_dir] - C1.iDimSize[X_dir];
                                }
                                if(P1.iPos[Y_dir] < 0){
                                        P1.dPos[Y_dir] = P1.dPos[Y_dir] + C1.dDimLength[Y_dir];
                                        P1.iPos[Y_dir] = P1.iPos[Y_dir] + C1.iDimSize[Y_dir];
                                }else if(P1.iPos[Y_dir] >= C1.iDimSize[Y_dir]){
                                        P1.dPos[Y_dir] = P1.dPos[Y_dir] - C1.dDimLength[Y_dir];
                                        P1.iPos[Y_dir] = P1.iPos[Y_dir] - C1.iDimSize[Y_dir];
                                }
                        }else{
                                if (  P1.iPos[X_dir] <  0 ||  P1.iPos[X_dir] >= C1.iDimSize[X_dir] ){
                                        P1.ParticleType = 0;
                                        continue;
                                }
                                if (   P1.iPos[Y_dir] < 0 ||  P1.iPos[Y_dir] >= C1.iDimSize[Y_dir] ){
                                        P1.ParticleType = 0;
                                        continue;
                                }
                        }

                        if(P1.iPos[Z_dir] < 0 || P1.iPos[Z_dir] >= C1.iDimSize[Z_dir]){
                                P1.ParticleType = 0;
                                continue;
                        }else{
                                itag =  P1.iPos[X_dir] + ( P1.iPos[Y_dir] + P1.iPos[Z_dir]*C1.iDimSize[Y_dir] )*C1.iDimSize[X_dir];
                        }


                        //--Calculate particle neighboring six points (particle size edge)
                        for ( int i = 0; i < 6 ; i++){
                                dPos_six_point[i][X_dir]   = P1.dPos[X_dir] + six_point_particle[i][X_dir]*dx_particle;
                                dPos_six_point[i][Y_dir]   = P1.dPos[Y_dir] + six_point_particle[i][Y_dir]*dy_particle;
                                dPos_six_point[i][Z_dir]   = P1.dPos[Z_dir] + six_point_particle[i][Z_dir]*dz_particle;
                                iPos_six_point[i][X_dir]     = int(   floor(dPos_six_point[i][X_dir]/C1.dDimLength_per_cell[X_dir] )  );
                                iPos_six_point[i][Y_dir]     = int(   floor( dPos_six_point[i][Y_dir]/C1.dDimLength_per_cell[Y_dir] ) );
                                iPos_six_point[i][Z_dir]     = int(   floor( dPos_six_point[i][Z_dir]/C1.dDimLength_per_cell[Z_dir] ) );
                                if(iPos_six_point[i][X_dir] < 0){
                                        dPos_six_point[i][X_dir] = dPos_six_point[i][X_dir] + C1.dDimLength[X_dir];
                                        iPos_six_point[i][X_dir] = iPos_six_point[i][X_dir] + C1.iDimSize[X_dir];
                                }else if(  iPos_six_point[i][X_dir] >= C1.iDimSize[X_dir]){
                                        dPos_six_point[i][X_dir] = dPos_six_point[i][X_dir] - C1.dDimLength[X_dir];
                                        iPos_six_point[i][X_dir] = iPos_six_point[i][X_dir] - C1.iDimSize[X_dir];
                                }
                                if( iPos_six_point[i][Y_dir] < 0){
                                        dPos_six_point[i][Y_dir] = dPos_six_point[i][Y_dir] + C1.dDimLength[Y_dir];
                                        iPos_six_point[i][Y_dir] = iPos_six_point[i][Y_dir] + C1.iDimSize[Y_dir];
                                }else if( iPos_six_point[i][Y_dir] >= C1.iDimSize[Y_dir]){
                                        dPos_six_point[i][Y_dir] = dPos_six_point[i][Y_dir] - C1.dDimLength[Y_dir];
                                        iPos_six_point[i][Y_dir] = iPos_six_point[i][Y_dir] - C1.iDimSize[Y_dir];
                                }
                                if(iPos_six_point[i][Z_dir] < 0 || iPos_six_point[i][Z_dir] >= C1.iDimSize[Z_dir]){
                                        P1.ParticleType = 0;
                                        continue;
                                }else{
                                        itag_six_point[i] =  iPos_six_point[i][X_dir] + ( iPos_six_point[i][Y_dir] + iPos_six_point[i][Z_dir]*C1.iDimSize[Y_dir] )*C1.iDimSize[X_dir];
                                }
                        }


                        //--count how many point for a seven-point molecule is on solid
                        count_point_on_solid =0;
                        for (int i = 0; i < 6 ; i++){
                                if (C1.iStatus[itag_six_point[i]] == iSubstrateStat || C1.iStatus[itag_six_point[i]] == iMaskStat){
                                        count_point_on_solid++;
                                }
                        }
                        if (C1.iStatus[itag] == iSubstrateStat || C1.iStatus[itag] == iMaskStat ){
                                count_point_on_solid++;
                        }



                        //--check if a particle collide on solid cell
                        if (count_point_on_solid > 1){
                                P1.dPos[X_dir] = old_dPos[X_dir];
                                P1.dPos[Y_dir] = old_dPos[Y_dir];
                                P1.dPos[Z_dir] = old_dPos[Z_dir];
                                P1.time_interval = 0.5 * P1.time_interval;
                        }else if ( count_point_on_solid == 1){
                                for (int i = 0; i < 6 ; i++){
                                        if (C1.iStatus[itag_six_point[i]] == iSubstrateStat || C1.iStatus[itag_six_point[i]] == iMaskStat){
                                                itag = itag_six_point[i];
                                        }
                                }

                                C1.surface_normal(searching_index, searching_number, itag, P1.iPos, P1.Vel,
                                                                          norm_surface_N, norm_reflected_V, &grazing_angle,  &incident_angle );

                                if (C1.iStatus[itag] == iMaskStat){

                                        if (P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType ){
                                                C1.IonMaskReaction(phys_sputter_prob, itag, P1.energy*Joule_to_eV, incident_angle, &ReactionExecution  );
                                        }

                                        P1.reflected_velocity_with_new_energy(norm_reflected_V,  &grazing_angle, P1.Vel);

                                        if ( P1.speed == 0){P1.ParticleType = 0;
                                        }else{                            P1.time_interval = dx/P1.speed;}


                                }else if ( C1.iStatus[itag] == iSubstrateStat ){

                                        if( P1.ParticleType == iClRadicalType ){
                                                C1.ClRadicalReaction(p0_ClRadicalReaction, itag, iNumMaterial, &ReactionExecution, &reaction_index);

                                                if (ReactionExecution == 0){
                                                        P1.speed = P1.setInitialSpeed( Temperature,   MassChlorine,   speed_cutoff_for_thermal_paricle[iClRadicalType]  );
                                                        P1.theta = P1.setReemitTheta(n_for_cosine_law); //--theta is with respect to surface normal
                                                        P1.phi = unif(generator)*2*PI;
                                                        if( norm_EmitParticle_V[X_dir] == 0 && norm_EmitParticle_V[Y_dir] == 0 ){
                                                                alpha = 0;
                                                        }else{
                                                                alpha = acos(norm_EmitParticle_V[X_dir] / sqrt( pow(norm_EmitParticle_V[X_dir],2.0)+pow(norm_EmitParticle_V[Y_dir],2.0)  )   );
                                                        }
                                                        beta = acos(norm_EmitParticle_V[Z_dir]);
                                                        P1.Vel[X_dir] = P1.speed*( cos(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+cos(beta)*sin(alpha)*cos(P1.theta)-sin(beta)*sin(P1.theta)*sin(P1.phi));
                                                        P1.Vel[Y_dir] = P1.speed*( sin(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)+cos(beta)*sin(P1.theta)*sin(P1.phi));
                                                        P1.Vel[Z_dir] = P1.speed*( -sin(alpha)*sin(P1.theta)*cos(P1.phi)+cos(alpha)*cos(P1.theta));
                                                        P1.time_interval = dx/P1.speed;
                                                }else{
                                                        P1.ParticleType = 0;
                                                }

                                        }else if ( P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType){

                                                if (P1.ParticleType == iClIonType){
                                                        C1.ClIonReaction(Eth_ClIonReaction, &E0_ClIonReaction, p0_ClIonReaction, type_ClIonReaction,
                                                                                               phys_sputter_prob, chem_sputter_prob, itag, iNumMaterial, P1.energy*Joule_to_eV,
                                                                                               incident_angle, &ReactionExecution, &ReflectedParticle, &EmitParticle, &reaction_index);
                                                }else if (P1.ParticleType == iCl2IonType){
                                                        C1.Cl2IonReaction(Eth_Cl2IonReaction, &E0_Cl2IonReaction, p0_Cl2IonReaction, type_Cl2IonReaction,
                                                                                                 phys_sputter_prob, chem_sputter_prob, itag, iNumMaterial, P1.energy*Joule_to_eV,
                                                                                                incident_angle, &ReactionExecution, &ReflectedParticle, &EmitParticle, &reaction_index);
                                                }else if (P1.ParticleType == iArIonType){
                                                        C1.ArIonReaction( Eth_ArIonReaction, &E0_ArIonReaction, p0_ArIonReaction, type_ArIonReaction,
                                                                                                phys_sputter_prob, chem_sputter_prob, itag, iNumMaterial, P1.energy*Joule_to_eV,
                                                                                                incident_angle, &ReactionExecution, &ReflectedParticle, &EmitParticle, &reaction_index);
                                                }

                                                if (ReactionExecution == 1 && EmitParticle != 0){
                                                        dPos_EmitParticle[X_dir] = P1.dPos[X_dir];
                                                        dPos_EmitParticle[Y_dir] = P1.dPos[Y_dir];
                                                        dPos_EmitParticle[Z_dir] = P1.dPos[Z_dir];
                                                        iPos_EmitParticle[X_dir] = P1.iPos[X_dir];
                                                        iPos_EmitParticle[Y_dir] = P1.iPos[Y_dir];
                                                        iPos_EmitParticle[Z_dir] = P1.iPos[Z_dir];
                                                        norm_EmitParticle_V[X_dir] = norm_surface_N[X_dir];
                                                        norm_EmitParticle_V[Y_dir] = norm_surface_N[Y_dir];
                                                        norm_EmitParticle_V[Z_dir] = norm_surface_N[Z_dir];
                                                }

                                                P1.ParticleType = ReflectedParticle;
                                                P1.reflected_velocity_with_new_energy(norm_reflected_V,  &grazing_angle, P1.Vel);

                                                if ( P1.speed == 0){  P1.ParticleType = 0;
                                                }else{                               P1.time_interval = dx/P1.speed;}

                                        }else if ( P1.ParticleType == iSiClgType || P1.ParticleType == iSiCl2gType || P1.ParticleType == iSiCl3gType){
                                                C1.redeposition(p0_redeposition, EmitParticle, itag, &ReactionExecution, &reaction_index);
                                                P1.ParticleType = 0;
                                        }
                                }
                        }
                }//--End of particle propagation loop
        }//--End of particle generation loop
        cout << "Total number of files generated : " << file_index <<endl;
        cout << "Total Particle Number : " << particleNumber << endl;
        file_index++;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumSiClxs, iNumMaterial, PrintSi, PrintSiCl, PrintSiCl2, PrintSiCl3,
                                                               directory+output_file_name, append, file_index );

        /*
        ofstream out1( directory+"number_of_reactions"   );
         for (int i = 0; i < 26 ; i++){
                out1 << reaction_number[i] << endl;
        }
        cout << "Reaction numbers are written to " << directory+"number_of_reactions" << endl;
        ofstream out2(  directory + "number_of_particles"  );
        for (int i = 0; i < 500 ; i++){
                out2 << num_particle[i] << endl;
        }
        */

        cout << "Particle numbers are written to " << directory+"number_of_particles" << endl;
}//--End of Program
