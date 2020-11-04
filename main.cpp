#include <stdlib.h>
#include <iomanip>
#include <string>
#include "etching.h"
#include "cell.h"
#include "particle.h"
#include "json.hpp"
#include "rand.h"

#define N 100
#define numBlocks 1
#define threadsPerBlock 100

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
        string inputJsonFile, tmpString, geometry, meshfile;
        double dx, dy, dz;
        int Nx, Ny, Nz, iSubstrateThickZ, iMaskThickZ, iTrenchWidthX, iMaskWidthX, iNumMaterial, iNumMask;
        int searching_radius, searching_range, searching_number ;


        string output_type, append, vtkDataType, directory, output_file_name;
        bool PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3 ;
        int output_file_number;


        vector<double> phys_sputter_prob, chem_sputter_prob;
        vector<double> p0_ClRadicalReaction, p0_redeposition, p0_ClIonReaction, p0_Cl2IonReaction, p0_ArIonReaction;
        vector<double> Eth_ClIonReaction, Eth_Cl2IonReaction, Eth_ArIonReaction;
        double E0_ClIonReaction, E0_Cl2IonReaction, E0_ArIonReaction;
        vector<int> type_ClIonReaction, type_Cl2IonReaction, type_ArIonReaction ;
        string IEADF_data_directory, ion_energy_data_file, ion_angle_data_file;
        string cumulativeflux_ClIon_file, cumulativeflux_Cl2Ion_file,  cumulativeflux_ArIon_file ,  boundary_condition;
        vector<double> cumulativeflux_ArIon, cumulativeflux_ClIon, cumulativeflux_Cl2Ion, ion_energy, ion_angle;
        bool NEUTRAL_THETA_GAUSSIAN, ION_THETA_GAUSSIAN;
        double Temperature, dClRadicalFlux, dSigFlux, dSiClgFlux, dSiCl2gFlux, dSiCl3gFlux, dSiCl4gFlux, dClIonFlux, dCl2IonFlux, dArIonType ;
        double neutral_gaussian_sigma, ion_gaussian_sigma, particle_size_factor, n_for_cosine_law;
        int timestep_number;



        double Lx, Ly, Lz ;
        int reaction_number[26] = {0};  //--for reaction index from 1 - 26







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

                //--mesh tag reading
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
	                    //--for structured grid, dx = dy = dz
	                    dx                                 = config["cell_size"] ; //--unit : m, x length per cell size
	                    dy                                 = config["cell_size"] ; //--unit : m, y length per cell size
	                    dz                                 = config["cell_size"] ; //--unit : m, z length per cell size
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
                        searching_range = 2*searching_radius+1;
                        searching_number = pow(  searching_range, 3);

                }else{
	                    cout << "No mesh tag in input file! " << endl ;
	                    exit( 1 ) ;
                }


                //--output_setting tag reading
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
                        PRINT_SI = config["PRINT_SI"];
                        PRINT_SICl = config["PRINT_SICl"];
                        PRINT_SICl2 = config["PRINT_SICl2"];
                        PRINT_SICl3 = config["PRINT_SICl3"];
                        output_file_number = config["output_file_number"] ;
                }else{
	                    cout << "No output_setting tag in input file! " << endl ;
	                    exit( 1 ) ;
                }


                //--surface_reactions tag reading
                if ( configuration_json.find( "surface_reactions" ) != configuration_json.end() ){
                        json config = configuration_json[ "surface_reactions" ] ;
                        phys_sputter_prob        = config["phys_sputter_prob"].get<std::vector<double>>();
                        chem_sputter_prob      = config["chem_sputter_prob"].get<std::vector<double>>();
                        p0_ClRadicalReaction = config["Cl_radical_reaction"]["probability"].get<std::vector<double>>();
                        p0_redeposition             = config["redeposition"]["probability"].get<std::vector<double>>();
                        p0_ClIonReaction         = config["Cl_ion_reaction"]["probability"].get<std::vector<double>>();
                        p0_Cl2IonReaction      = config["Cl2_ion_reaction"]["probability"].get<std::vector<double>>();
                        p0_ArIonReaction        =  config["Ar_ion_reaction"]["probability"].get<std::vector<double>>();
                        Eth_ClIonReaction       =  config["Cl_ion_reaction"]["activation_energy"].get<std::vector<double>>();
                        Eth_Cl2IonReaction    =  config["Cl2_ion_reaction"]["activation_energy"].get<std::vector<double>>();
                        Eth_ArIonReaction      =  config["Ar_ion_reaction"]["activation_energy"].get<std::vector<double>>();
                        E0_ClIonReaction        =  config["Cl_ion_reaction"]["reference_energy"];
                        E0_Cl2IonReaction     = config["Cl2_ion_reaction"]["reference_energy"];
                        E0_ArIonReaction       = config["Ar_ion_reaction"]["reference_energy"];
                        type_ClIonReaction    = config["Cl_ion_reaction"]["reaction_type"].get<std::vector<int>>();
                        type_Cl2IonReaction = config["Cl2_ion_reaction"]["reaction_type"].get<std::vector<int>>();
                        type_ArIonReaction   = config["Ar_ion_reaction"]["reaction_type"].get<std::vector<int>>();
                }else{
                        cout << "No surface_reactions tag in input file! " << endl ;
	                    exit( 1 ) ;
                }



                //--simulation_conditions tag reading
                if ( configuration_json.find( "simulation_conditions" ) != configuration_json.end() ){
                        json config = configuration_json[ "simulation_conditions" ] ;
                        //cout << config << endl;
                        if ( config.find( "boundary_condition" ) != config.end() ){
                                boundary_condition = config["boundary_condition"];
	                    }else{
	                            boundary_condition = "periodic" ;
	                    }

                        IEADF_data_directory = config["IEADF_data_directory"];
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
	                    NEUTRAL_THETA_GAUSSIAN = config["NEUTRAL_THETA_GAUSSIAN"];
	                    ION_THETA_GAUSSIAN = config["ION_THETA_GAUSSIAN"];
	                    Temperature = config["Temperature"];
	                    dClRadicalFlux =config["dClRadicalFlux"];
	                    dSigFlux = config["dSigFlux"];
	                    dSiClgFlux = config["dSiClgFlux"];
	                    dSiCl2gFlux = config["dSiCl2gFlux"];
	                    dSiCl3gFlux = config["dSiCl3gFlux"];
	                    dSiCl4gFlux = config["dSiCl4gFlux"];
                        dClIonFlux = config["dClIonFlux"];
                        dCl2IonFlux = config["dCl2IonFlux"];
                        dArIonType = config["dArIonFlux"];


                        neutral_gaussian_sigma = config["neutral_gaussian_sigma"];
                        ion_gaussian_sigma = config["ion_gaussian_sigma"];
                        particle_size_factor = config["particle_size_factor"];

                        n_for_cosine_law = config["n_for_cosine_law"];
                        if ( config.find( "timestep_number" ) != config.end() ){
	                            timestep_number = config[ "timestep_number" ] ;
                        }
                }else{
                        cout << "No simulation_conditions tag in input file! " << endl ;
	                    exit( 1 ) ;
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
                cout << "Output setting:" << endl ;
                cout << "  form                         =    " << output_type << endl ;
                cout << "  vtkDataType                  =    " << vtkDataType << endl ;
                cout << "  output_file_directory        =    "<< directory << endl ;
                cout << "  output_file_name             =    "<< output_file_name << endl ;
                cout << "  output_file_number           =    "<< output_file_number << endl;
                cout << "  PRINT_SI                     =    " << PRINT_SI << endl ;
                cout << "  PRINT_SICl                   =    " << PRINT_SICl << endl ;
                cout << "  PRINT_SICl2                  =    " << PRINT_SICl2 << endl ;
                cout << "  PRINT_SICl3                  =    " << PRINT_SICl3 << endl ;
                cout << endl;
                cout << "Simulation conditions:" << endl ;
                cout << "  timestep_number              =    " << timestep_number << endl ;
                cout << "  Boundary condition           =    " << boundary_condition << endl;
                cout << "  Temperature                  =    "  << Temperature << endl;
                cout << "  Cl radical flux              =    " << dClRadicalFlux << endl;
	            cout << "  Cl ion flux                  =    " <<        dClIonFlux << endl;
                cout << "  Cl2+ ion flux                =    " << dCl2IonFlux << endl;
                cout <<  "  Ar ion flux                  =    " <<    dArIonType << endl;
                cout << "  IEADF_data_directory         =    " << IEADF_data_directory << endl;
                cout << "  NEUTRAL_THETA_GAUSSIAN       =    " << NEUTRAL_THETA_GAUSSIAN << endl;
                cout << "  neutral_gaussian_sigma       =    " << neutral_gaussian_sigma << endl;
                cout << "  ION_THETA_GAUSSIAN           =    " << ION_THETA_GAUSSIAN << endl;
                cout << "  ion_gaussian_sigma           =    " << ion_gaussian_sigma << endl;
                cout << "  n for cosine law             =    " <<  n_for_cosine_law <<endl;
                cout << "  particle size factor         =    " <<  particle_size_factor <<endl;
                cout << endl;
        }//--End of reading and printing json file















        /* Initialization for cell */
        cell C1;
        if ( meshfile == "none" ){
                //dy = dx; //--y length per cell size
                //dz = dx; //--z length per cell size
                Lx = Nx*dx;
                Ly = Ny*dy;
                Lz = Nz*dz;
                //dSubstrateThickZ = iSubstrateThickZ*dz;
                //dMaskThickZ = iMaskThickZ*dz;
                //dTrenchWidthX = iTrenchWidthX*dx;
                //dMaskWidthX = iMaskWidthX*dx;
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
                                string meshDATASET =  buffer_string_line;
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




        /*random number generation*/
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution


        /* Pre-calculation of speed cutoffe */
        int ThermalParticleTypes = 6; //--Total number of thermal species, here we include Cl radical, Si(g), SiCl(g), SiCl2(g), and SiCl3(g)
        double speed_cutoff_for_thermal_paricle [ThermalParticleTypes]; //--Array to store pre-calculated cutoff speed
        speed_cutoff_for_thermal_paricle[iClRadicalType] = calc_v_cut(  MassChlorine, Temperature  );
        speed_cutoff_for_thermal_paricle[iSigType]              = calc_v_cut(  MassSilicon, Temperature  );
        speed_cutoff_for_thermal_paricle[iSiClgType]          = calc_v_cut(  MassSilicon+MassChlorine, Temperature  );
        speed_cutoff_for_thermal_paricle[iSiCl2gType]       = calc_v_cut(  MassSilicon+2*MassChlorine, Temperature  );
        speed_cutoff_for_thermal_paricle[iSiCl3gType]       = calc_v_cut(  MassSilicon+3*MassChlorine, Temperature  );


         /* Pre-calculation for particle generation probability */
        double TotalFlux = dClRadicalFlux+dSigFlux+dSiClgFlux+dSiCl2gFlux+dSiCl3gFlux+dSiCl4gFlux+dClIonFlux+dCl2IonFlux+dArIonType;
        vector<double> ParticleProb_for_incident_particle;
        ParticleProb_for_incident_particle.push_back(dClRadicalFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dSigFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dSiClgFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dSiCl2gFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dSiCl3gFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dSiCl4gFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dClIonFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dCl2IonFlux/TotalFlux);
        ParticleProb_for_incident_particle.push_back(dArIonType/TotalFlux);


        /*Pre-calculation of particle size and initialization of six point particle*/
        double six_point_particle [6][3] = {  {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1} };
        double dx_particle, dy_particle, dz_particle;
        dx_particle = particle_size_factor*dx;
        dy_particle = particle_size_factor*dy;
        dz_particle = particle_size_factor*dz;


        /*Pre-calculation for surface sites*/
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


        /*Pre-calculation of total particle number, time interval,  and output frequency*/
        double Gamma_t = (dClRadicalFlux + dClIonFlux + dCl2IonFlux + dArIonType)*1E4 ; //--convert unit from cm^-2 to m^-2
        double flux_area = Lx*Ly;
        double total_time;
        double time_interval = 1/Gamma_t/flux_area;
        int TotalParticle;
        string yes_or_no = "yes";
        cout  << "The area of incidence (A) = " << flux_area << " m^2" << endl;
        cout << "The total flux (Gamma_t) = "<< Gamma_t << " m^-2 s^-1" << endl;
        //cout << "The number of atoms in a cell (Ns) = " << iNumMaterial << endl;
        cout << "The time interval between steps (dt) = " << "1/(Gamma_t*A) = " << time_interval << endl;
        do{
                cout << "How many seconds in real time would you like to simulate? " ;
                cin >> total_time;
                TotalParticle = int(total_time/time_interval);
                cout << "The total number of steps in the simulation will be : " << TotalParticle << endl;
                cout << "Is that OK with you (yes/no)? " ;
                cin >> yes_or_no;
        } while ( yes_or_no == "no");
        int file_index = 0;
        int frequency = int(TotalParticle/output_file_number);
        int particleNumber = 0;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                     PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+output_file_name, append, file_index );




        /*
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
        */

        #pragma omp parallel for
        for (int indexParticle = 0;   indexParticle < TotalParticle;   indexParticle++){



                particle P1;                                            //--Generate a particle
                P1.setInitialType(ParticleProb_for_incident_particle);  //--choose a particular type of particle
                double rdd1, rdd2, rdd3;                                //--random number for selecting ion energy and angle
                int ReactionExecution = 0;                              //--determine if a reaction happen
                int EmitParticle = 0;                                   //--index for emitted particle such as SiClx(g)
                int EmitTimes = 0;
                int ReflectedParticle = 0;                              //--index for original reflected particle such as Ar*, Cl*
                int reaction_index = 0;                                 //--index for a particular reaction formula
                double norm_surface_N [3];                              //--normalized surface normal vector (x, y, z)
                double norm_reflected_V [3];                            //--normalized reflected velocity vector (x, y, z)
                double reflected_velocity [3];                          //--reflected velocity (x, y, z)
                double grazing_angle = 0;                               //--angle between surface and velocity
                double incident_angle = 0;                              //--angle between normal and velocity
                int old_iPos [3];                                       //--to record the iPos (x, y, z) before propogation
                double old_dPos [3];                                    //--to record the dPos (x, y, z) before propagation
                int itag;                                               //--itag of particle center
                int old_itag;
                int itag_six_point [6] ;                                //--itag of particle's six neighboring points
                double dPos_six_point [6][3] ;                          //--dPos (x, y, z) of particle's six neighboring points
                int iPos_six_point [6][3] ;                             //--iPos (x, y, z) of particle's siz neighboring points
                int count_point_on_solid;                               //--to count how many point of a seven-point molecule is on solid cell
                double norm_EmitParticle_V [3];                         //--normalized velocity for emitted particle
                double dPos_EmitParticle [3];                           //--dPos for emitted particle
                int iPos_EmitParticle [3];                              //--iPos for emitted particle
                double alpha;                                           //--azimusal angle of surface normal vector
                double beta;                                            //--polar angle of surface normal vector;

                //--Count total number of particle and write to a file in a vtk format
                #pragma omp critical
                {
                        particleNumber++;
                        if(  particleNumber%frequency == 0  ){
                                file_index++;
                                if (file_index == output_file_number){
                                        PRINT_SI = true;
                                        PRINT_SICl = true;
                                        PRINT_SICl2 = true;
                                        PRINT_SICl3 = true;
                                }
                                write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                                                PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+output_file_name, append, file_index );
                        }
                }

                if (P1.ParticleType == 0)    continue;

                if ( P1.ParticleType == iClRadicalType ){
                        P1.mass = MassChlorine;
                        P1.speed = P1.setInitialSpeed( Temperature,   MassChlorine,   speed_cutoff_for_thermal_paricle[iClRadicalType]  );
                        P1.energy = 0.5 * P1.mass * P1.speed * P1.speed; //--unit: Joule
                        if ( NEUTRAL_THETA_GAUSSIAN == true){
                                P1.theta = P1.setInitialTheta_by_Gaussian(neutral_gaussian_sigma); //--unit: radian
                        }else if ( NEUTRAL_THETA_GAUSSIAN == false){
                                P1.theta = unif(generator)*PI/2+PI/2;          //--generate the random number between PI/2 to PI for theta
                        }
                }else if ( P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType ){
                        if ( P1.ParticleType == iClIonType){
                                P1.mass = MassChlorine;
                        }else if ( P1.ParticleType == iCl2IonType){
                                P1.mass = MassChlorine*2;
                        }else if (  P1.ParticleType == iArIonType){
                                P1.mass = MassArgon;
                        }
                        rdd1 = unif(generator);
                        rdd2 = unif(generator);
                        rdd3 = unif(generator);
                        if ( ION_THETA_GAUSSIAN == true){
                                P1.theta = P1.setInitialTheta_by_Gaussian(ion_gaussian_sigma); //--unit: radian
                        }else if ( ION_THETA_GAUSSIAN == false){
                                for(int i = 0;    i < cumulativeflux_ClIon.size()-1;     i++){
                                        if  (  rdd1 >= cumulativeflux_ClIon[i] && rdd1 < cumulativeflux_ClIon[i+1]  ){
                                                P1.theta =(  ion_angle[i] + rdd2*(ion_angle[i+1] - ion_angle[i])  )*PI/180+PI; //--unit: radian
                                                break;
                                        }
                                }
                        }
                        for(int i = 0;    i < cumulativeflux_ClIon.size()-1;     i++){
                                if  (  rdd1 >= cumulativeflux_ClIon[i] && rdd1 < cumulativeflux_ClIon[i+1]  ){
                                        P1.energy = (  ion_energy[i] + rdd3*(ion_energy[i+1] - ion_energy[i]) )/Joule_to_eV;//--unit: Joule
                                        break;
                                }
                        }
                        P1.speed = sqrt(2*P1.energy/P1.mass);

                }
                P1.phi = unif(generator)*2*PI; //--unit: radian
                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                P1.time_interval = dx/P1.speed;
                P1.Vel[X_dir] = P1.speed*sin(P1.theta)*cos(P1.phi);
                P1.Vel[Y_dir] = P1.speed*sin(P1.theta)*sin(P1.phi);
                P1.Vel[Z_dir] = P1.speed*cos(P1.theta);


                for(int indexTimeStep = 0; indexTimeStep < timestep_number; indexTimeStep++){

                        //--check if particle has been deactivated
                        if(P1.ParticleType == 0){

                                if (EmitParticle == 0){
                                        break;
                                }else{
                                        //if (EmitTimes > 0)    break;
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
                                                //alpha = acos(norm_EmitParticle_V[X_dir] / sqrt( pow(norm_EmitParticle_V[X_dir],2.0)+pow(norm_EmitParticle_V[Y_dir],2.0)  )   );

                                                if ( norm_EmitParticle_V[X_dir] >= 0){
                                                        alpha = acos(-norm_EmitParticle_V[Y_dir] / sqrt( pow(norm_EmitParticle_V[X_dir],2.0)+pow(norm_EmitParticle_V[Y_dir],2.0)  )   );
                                                }else if ( norm_EmitParticle_V[X_dir] < 0){
                                                        alpha = -acos(-norm_EmitParticle_V[Y_dir] / sqrt( pow(norm_EmitParticle_V[X_dir],2.0)+pow(norm_EmitParticle_V[Y_dir],2.0)  )   );
                                                }

                                        }
                                        beta = acos(norm_EmitParticle_V[Z_dir]);

                                        //--Use the formula that I derived
                                        P1.Vel[X_dir] = P1.speed*( -sin(alpha)*cos(beta)*sin(P1.theta)*sin(P1.phi)+cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(alpha)*sin(beta)*cos(P1.theta)                             );
                                        P1.Vel[Y_dir] = P1.speed*(   cos(alpha)*cos(beta)*sin(P1.theta)*sin(P1.phi)+sin(alpha)*sin(P1.theta)*cos(P1.phi)-cos(alpha)*sin(beta)*cos(P1.theta)                              );
                                        P1.Vel[Z_dir] = P1.speed*(   sin(beta)*sin(P1.theta)*sin(P1.phi)+cos(beta)*cos(P1.theta)                                                    );


                                        //--Use the formula proposed by Kushner in
                                        /*
                                        P1.Vel[X_dir] = P1.speed*(   cos(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+cos(beta)*sin(alpha)*cos(P1.theta)-sin(beta)*sin(P1.theta)*sin(P1.phi)                             );
                                        P1.Vel[Y_dir] = P1.speed*(   sin(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)+cos(beta)*sin(P1.theta)*sin(P1.phi)                              );
                                        P1.Vel[Z_dir] = P1.speed*( -sin(alpha)*sin(P1.theta)*cos(P1.phi)+cos(alpha)*cos(P1.theta)                                                    );
                                        */


                                        P1.time_interval = dx/P1.speed;
                                        EmitParticle = 0;
                                        //EmitTimes++;
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

                                C1.surface_normal(searching_index, searching_number, itag, P1.iPos, P1.Vel, norm_surface_N, norm_reflected_V, &grazing_angle,  &incident_angle );

                                if (C1.iStatus[itag] == iMaskStat){

                                        if (P1.ParticleType == iClRadicalType){
                                                P1.Vel[X_dir] = P1.speed*norm_reflected_V[X_dir];
                                                P1.Vel[Y_dir] = P1.speed*norm_reflected_V[Y_dir];
                                                P1.Vel[Z_dir] = P1.speed*norm_reflected_V[Z_dir];
                                        }


                                        if (P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType ){
                                                break; //--for test
                                                C1.IonMaskReaction(phys_sputter_prob, itag, P1.energy*Joule_to_eV, incident_angle, &ReactionExecution  );
                                                P1.reflected_velocity_with_new_energy(norm_reflected_V,  &grazing_angle, P1.Vel);
                                        }


                                        if (P1.ParticleType == iSiClgType || P1.ParticleType == iSiCl2gType || P1.ParticleType == iSiCl3gType ){
                                                break; //--for test
                                        }


                                        if ( P1.speed == 0){
                                                P1.ParticleType = 0;
                                        }else{
                                                P1.time_interval = dx/P1.speed;
                                        }

                                        continue;


                                }else if ( C1.iStatus[itag] == iSubstrateStat ){

                                        if( P1.ParticleType == iClRadicalType ){
                                                C1.ClRadicalReaction(p0_ClRadicalReaction, itag, iNumMaterial, &ReactionExecution, &reaction_index);

                                                if (ReactionExecution == 0){
                                                        /*  for test
                                                        cout << "Neutral radicals hit the substrate but no reaction happens !" ;
                                                        cin.get();
                                                        */
                                                        P1.speed = P1.setInitialSpeed( Temperature,   MassChlorine,   speed_cutoff_for_thermal_paricle[iClRadicalType]  );
                                                        P1.theta = P1.setReemitTheta(n_for_cosine_law); //--theta is with respect to surface normal
                                                        P1.phi = unif(generator)*2*PI;
                                                        if( norm_surface_N[X_dir] == 0 && norm_surface_N[Y_dir] == 0 ){
                                                                alpha = 0;
                                                        }else{
                                                                 //alpha = acos(norm_surface_N[X_dir] / sqrt( pow(norm_surface_N[X_dir],2.0)+pow(norm_surface_N[Y_dir],2.0)  )   );


                                                                if ( norm_surface_N[X_dir] >= 0){
                                                                        alpha = acos(-norm_surface_N[Y_dir] / sqrt( pow(norm_surface_N[X_dir],2.0)+pow(norm_surface_N[Y_dir],2.0)  )   );
                                                                }else if ( norm_surface_N[X_dir] < 0){
                                                                        alpha = -acos(-norm_surface_N[Y_dir] / sqrt( pow(norm_surface_N[X_dir],2.0)+pow(norm_surface_N[Y_dir],2.0)  )   );
                                                                }

                                                        }
                                                        beta = acos(norm_surface_N[Z_dir]);

                                                        //--Use the formula proposed by Kushner in
                                                        /*
                                                        P1.Vel[X_dir] = P1.speed*(   cos(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+cos(beta)*sin(alpha)*cos(P1.theta)-sin(beta)*sin(P1.theta)*sin(P1.phi)                             );
                                                        P1.Vel[Y_dir] = P1.speed*(   sin(beta)*cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)+cos(beta)*sin(P1.theta)*sin(P1.phi)                              );
                                                        P1.Vel[Z_dir] = P1.speed*( -sin(alpha)*sin(P1.theta)*cos(P1.phi)+cos(alpha)*cos(P1.theta)                                                    );
                                                         */

                                                        //--Use the formula that I derived

                                                        P1.Vel[X_dir] = P1.speed*( -cos(beta)*sin(alpha)*sin(P1.theta)*sin(P1.phi)+cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)                             );
                                                        P1.Vel[Y_dir] = P1.speed*( cos(beta)*cos(alpha)*sin(P1.theta)*sin(P1.phi)+sin(alpha)*sin(P1.theta)*cos(P1.phi)-sin(beta)*cos(alpha)*cos(P1.theta)                              );
                                                        P1.Vel[Z_dir] = P1.speed*( sin(beta)*sin(P1.theta)*sin(P1.phi)+cos(beta)*cos(P1.theta)    );

                                                        P1.time_interval = dx/P1.speed;
                                                }else if (ReactionExecution == 1){
                                                        /*
                                                        cout << "Neutral radicals hit the substrate and reaction happens !" ;
                                                        cin.get();
                                                        */
                                                        P1.ParticleType = 0;
                                                }
                                                continue;
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
                                                        /*
                                                        cout << "Particle with energetic type hit the substrate and reaction happens with emitted particle !";
                                                        cin.get();
                                                        */
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

                                                if ( P1.speed == 0){
                                                        P1.ParticleType = 0;
                                                }else{
                                                        P1.time_interval = dx/P1.speed;
                                                }

                                                continue;

                                        }else if ( P1.ParticleType == iSiClgType || P1.ParticleType == iSiCl2gType || P1.ParticleType == iSiCl3gType){
                                                C1.redeposition(p0_redeposition, P1.ParticleType, itag, &ReactionExecution, &reaction_index);

                                                if (ReactionExecution == 0){
                                                        /*
                                                        cout << "Particle with thermal type hit the substrate but no reaction happens ! " << endl;
                                                        cin.get();
                                                        */
                                                        if ( P1.ParticleType == iSiClgType){
                                                                P1.speed = P1.setInitialSpeed( Temperature, MassSilicon+MassChlorine, speed_cutoff_for_thermal_paricle[iSiCl2gType]  );
                                                        }else if ( P1.ParticleType == iSiCl2gType){
                                                                P1.speed = P1.setInitialSpeed( Temperature, MassSilicon+2*MassChlorine, speed_cutoff_for_thermal_paricle[iSiCl2gType]  );
                                                        }else if ( P1.ParticleType == iSiCl3gType){
                                                                P1.speed = P1.setInitialSpeed( Temperature, MassSilicon+3*MassChlorine, speed_cutoff_for_thermal_paricle[iSiCl3gType]  );
                                                        }

                                                        P1.theta = P1.setReemitTheta(n_for_cosine_law); //--theta is with respect to surface normal
                                                        P1.phi = unif(generator)*2*PI;
                                                        if( norm_surface_N[X_dir] == 0 && norm_surface_N[Y_dir] == 0 ){
                                                                alpha = 0;
                                                        }else{
                                                                if ( norm_surface_N[X_dir] >= 0){
                                                                        alpha = acos(-norm_surface_N[Y_dir] / sqrt( pow(norm_surface_N[X_dir],2.0)+pow(norm_surface_N[Y_dir],2.0)  )   );
                                                                }else if ( norm_surface_N[X_dir] < 0){
                                                                        alpha = -acos(-norm_surface_N[Y_dir] / sqrt( pow(norm_surface_N[X_dir],2.0)+pow(norm_surface_N[Y_dir],2.0)  )   );
                                                                }

                                                        }
                                                        beta = acos(norm_surface_N[Z_dir]);

                                                        //--Use the formula that I derived
                                                        P1.Vel[X_dir] = P1.speed*( -cos(beta)*sin(alpha)*sin(P1.theta)*sin(P1.phi)+cos(alpha)*sin(P1.theta)*cos(P1.phi)+sin(beta)*sin(alpha)*cos(P1.theta)                             );
                                                        P1.Vel[Y_dir] = P1.speed*( cos(beta)*cos(alpha)*sin(P1.theta)*sin(P1.phi)+sin(alpha)*sin(P1.theta)*cos(P1.phi)-sin(beta)*cos(alpha)*cos(P1.theta)                              );
                                                        P1.Vel[Z_dir] = P1.speed*( sin(beta)*sin(P1.theta)*sin(P1.phi)+cos(beta)*cos(P1.theta)    );

                                                        P1.time_interval = dx/P1.speed;
                                                }else if (ReactionExecution == 1){
                                                        P1.ParticleType = 0;
                                                }

                                                continue;
                                        }//--End of emittion particle ( if clause)
                                }//--End of substrate (if clause)
                        }//--End of count_point_on_solid (if clause)
                }//--End of particle propagation loop (for loop)
        }//--End of particle generation loop (for loop)
        cout << "Total number of files generated : " << file_index <<endl;
        cout << "Total Particle Number : " << particleNumber << endl;
        file_index++;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                     PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+output_file_name, append, file_index );

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
