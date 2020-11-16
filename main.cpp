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
        //--parameterf read from mesh
        string inputJsonFile, tmpString, geometry, meshfile;
        double dx, dy, dz;
        int Nx, Ny, Nz, iSubstrateThickZ, iMaskThickZ, iTrenchWidthX, iMaskWidthX, iMaskLength, iGapWidthX, iGapWidthY,
               iNumMaterial, iNumMask;
        int SurfaceSearchingRadius, SurfaceSearchingRange, SurfaceSearchingNumber ;

        //--parameterd read from output setting
        string output_type, append, vtkDataType, directory, OutputFilename;
        bool PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3 ;
        int OutputFileNumber;

        //--parameters read from simulation conditions
        string IEADFDataDirectory, IonEnergyDataFilename, IonAngleDataFilename;
        string cumulativefluxClIonFilename, cumulativefluxCl2IonFilename,  cumulativefluxArIonFilename, BoundaryCondition;
        vector<double>  cumulativeflux_ClIon, cumulativeflux_Cl2Ion, cumulativeflux_ArIon, IonEnergy, IonAngle;
        bool NEUTRAL_THETA_SCALING, ION_THETA_SCALING;
        double Temperature, dClRadicalFlux, dClIonFlux, dCl2IonFlux, dArIonType, NeutralThetaScalingFactor, IonThetaScalingFactor,
                        ParticleSizeFactor, ReemissionCosineLawPower, TotalRealTime;
        int PropagationTimestepNumber;

        //--parameters read from surface ractions
        vector<double> PhysSputterProb, ChemSputterProb;
        vector<double> p0_ClRadicalReaction, p0_redeposition, p0_ClIonReaction, p0_Cl2IonReaction, p0_ArIonReaction;
        vector<double> Eth_ClIonReaction, Eth_Cl2IonReaction, Eth_ArIonReaction;  //--Eth is the activation energy for ion-enhanced reactions
        double E0_ClIonReaction, E0_Cl2IonReaction, E0_ArIonReaction;  //--E0 is the reference energy
        vector<int> type_ClIonReaction, type_Cl2IonReaction, type_ArIonReaction ;  //--Type of reactions: either physical or chemical sputtering






        double Lx, Ly, Lz ;
        int reaction_number[26] = {0};  //--for reaction index from 1 - 26
        double epsilon_0 = 0.0;  //--unit: eV
        double epsilon_s = 50;  //--unit: eV
        double theta_0 = 30;  //--unit: degree
        double gamma_0 = 0.85;  //--a scaling factor, no unit






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
	                    dx                                 = config["cell_size"] ; //--unit : nm, x length per cell size
	                    dy                                 = config["cell_size"] ; //--unit : nm, y length per cell size
	                    dz                                 = config["cell_size"] ; //--unit : nm, z length per cell size
	                    dx *= 1E-9;                                                              //--convert unit from nm to m
	                    dy *= 1E-9;
	                    dz *= 1E-9;
	                    Nx                                = config["domain_size_x"]; //--number of cells in x direction
                        Ny                                = config["domain_size_y"]; //--number of cells in y direction
                        Nz                                = config["domain_size_z"]; //--number of cells in z direction
                        Lx = Nx*dx;
                        Ly = Ny*dy;
                        Lz = Nz*dz;
                        iSubstrateThickZ = config["substrate_thickness"];
                        iMaskThickZ           = config["mask_thickness"];
                        iTrenchWidthX      = config["trench_width"];
                        iMaskWidthX         = config["mask_width"];
                        iMaskLength          = config["mask_length"];
                        iGapWidthX           = config["gap_width_x"];
                        iGapWidthY           = config["gap_width_y"];
                        iNumMaterial        = config["atom_number_in_material"];
                        iNumMask              = config["atom_number_in_mask"];
                        SurfaceSearchingRadius = config["surface_searching_radius"] ;
                        SurfaceSearchingRange = 2*SurfaceSearchingRadius+1;
                        SurfaceSearchingNumber = pow(  SurfaceSearchingRange, 3);

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
	                    OutputFilename = config["output_file_name"];
                        OutputFileNumber = config["output_file_number"] ;
                        PRINT_SI = config["PRINT_SI"];
                        PRINT_SICl = config["PRINT_SICl"];
                        PRINT_SICl2 = config["PRINT_SICl2"];
                        PRINT_SICl3 = config["PRINT_SICl3"];

                }else{
	                    cout << "No output_setting tag in input file! " << endl ;
	                    exit( 1 ) ;
                }


                //--simulation_conditions tag reading
                if ( configuration_json.find( "simulation_conditions" ) != configuration_json.end() ){
                        json config = configuration_json[ "simulation_conditions" ] ;
                        //cout << config << endl;
                        if ( config.find( "boundary_condition" ) != config.end() ){
                                BoundaryCondition = config["boundary_condition"];
	                    }else{
	                            BoundaryCondition = "periodic" ;
	                    }

                        IEADFDataDirectory = config["IEADF_data_directory"];
                        IonEnergyDataFilename =            IEADFDataDirectory + "energy.dat";
                        IonAngleDataFilename =               IEADFDataDirectory + "angle.dat";
                        cumulativefluxClIonFilename =   IEADFDataDirectory + "cumulativefluxClIon_data.dat";
                        cumulativefluxCl2IonFilename = IEADFDataDirectory + "cumulativefluxCl2Ion_data.dat";
                        cumulativefluxArIonFilename =   IEADFDataDirectory + "cumulativefluxArIon_data.dat";
                        IonEnergy                         = read_data(   IonEnergyDataFilename    );
                        IonAngle                            = read_data(   IonAngleDataFilename    );
                        cumulativeflux_ClIon    = read_data(   cumulativefluxClIonFilename   );
                        cumulativeflux_Cl2Ion = read_data(   cumulativefluxCl2IonFilename  );
                        cumulativeflux_ArIon   = read_data(   cumulativefluxArIonFilename   );
	                    NEUTRAL_THETA_SCALING = config["NEUTRAL_THETA_SCALING"];
	                    ION_THETA_SCALING = config["ION_THETA_SCALING"];
	                    Temperature = config["Temperature"];
	                    dClRadicalFlux =config["dClRadicalFlux"];
                        dClIonFlux = config["dClIonFlux"];
                        dCl2IonFlux = config["dCl2IonFlux"];
                        dArIonType = config["dArIonFlux"];
                        NeutralThetaScalingFactor = config["neutral_theta_scaling_factor"];
                        IonThetaScalingFactor = config["ion_theta_scaling_factor"];
                        ParticleSizeFactor = config["particle_size_factor"];
                        ReemissionCosineLawPower = config["reemission_cosine_law_power"];
                        if ( config.find( "propagation_timestep_number" ) != config.end() ){
	                            PropagationTimestepNumber = config[ "propagation_timestep_number" ] ;
                        }
                }else{
                        cout << "No simulation_conditions tag in input file! " << endl ;
	                    exit( 1 ) ;
                }


                //--surface_reactions tag reading
                if ( configuration_json.find( "surface_reactions" ) != configuration_json.end() ){
                        json config = configuration_json[ "surface_reactions" ] ;
                        PhysSputterProb        = config["phys_sputter_prob"].get<std::vector<double>>();
                        ChemSputterProb      = config["chem_sputter_prob"].get<std::vector<double>>();
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
                cout << "  mask_length                      =    " << iMaskLength << endl;
                cout << "  gap_width_x                  =    " << iGapWidthX << endl;
                cout << "  gap_width_y                  =    " << iGapWidthY << endl;
                cout << "  atom_number_in_material      =    " << iNumMaterial << endl;
                cout << "  atom_number_in_mask          =    " << iNumMask << endl;
                cout << "  surface_searching_radius     =    " << SurfaceSearchingRadius << endl;
                cout << endl;
                cout << "Output setting:" << endl ;
                cout << "  form                         =    " << output_type << endl ;
                cout << "  vtkDataType                  =    " << vtkDataType << endl ;
                cout << "  output_file_directory        =    "<< directory << endl ;
                cout << "  output_file_name             =    "<< OutputFilename << endl ;
                cout << "  output_file_number           =    "<< OutputFileNumber << endl;
                cout << "  PRINT_SI                     =    " << PRINT_SI << endl ;
                cout << "  PRINT_SICl                   =    " << PRINT_SICl << endl ;
                cout << "  PRINT_SICl2                  =    " << PRINT_SICl2 << endl ;
                cout << "  PRINT_SICl3                  =    " << PRINT_SICl3 << endl ;
                cout << endl;
                cout << "Simulation conditions:" << endl ;
                cout << "  total_real_time   =    "  << TotalRealTime << endl;
                cout << "  propagation_timestep_number  =    " << PropagationTimestepNumber << endl ;
                cout << "  Boundary condition           =    " << BoundaryCondition << endl;
                cout << "  Temperature                  =    "  << Temperature << endl;
                cout << "  Cl radical flux              =    " << dClRadicalFlux << endl;
	            cout << "  Cl ion flux                  =    " <<        dClIonFlux << endl;
                cout << "  Cl2+ ion flux                =    " << dCl2IonFlux << endl;
                cout <<  "  Ar ion flux                  =    " <<    dArIonType << endl;
                cout << "  IEADF_data_directory         =    " << IEADFDataDirectory << endl;
                cout << "  NEUTRAL_THETA_SCALING        =    " << NEUTRAL_THETA_SCALING << endl;
                cout << "  ION_THETA_SCALING            =    " << ION_THETA_SCALING << endl;
                cout << "  neutral_theta_scaling_factor =    " << NeutralThetaScalingFactor << endl;
                cout << "  ion_theta_scaling_factor     =    " << IonThetaScalingFactor << endl;
                cout << "  reemission_cosine_law_power  =    " <<  ReemissionCosineLawPower <<endl;
                cout << "  particle_size_factor         =    " <<  ParticleSizeFactor <<endl;
                cout << endl;
        }//--End of reading and printing json file















        /* Initialization for cell */
        cell C1;
        if ( meshfile == "single_long_trench" ){
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
        }else if ( meshfile == "array_set" ){
                C1.initial(Nx, Ny, Nz, Lx, Ly, Lz);
                for(int iz = 0; iz < Nz; iz++){
                        for(int iy = 0; iy < Ny; iy++){
                                for(int ix = 0; ix < Nx; ix++){
                                        int itag =  ix + ( iy + iz*Ny )*Nx;
                                        if(  iz < iSubstrateThickZ ){
                                                //--setting for substrate
                                                C1.setStatus(itag, iSubstrateStat, iNumMaterial);     //--Si has 8 atoms per unit cell with the dimension of 0.54 nm
                                        }else if (  iz >= iSubstrateThickZ && iz < (iSubstrateThickZ + iMaskThickZ)  ) {
                                                if ( ix < iGapWidthX/2 ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2 && ix < iGapWidthX/2+iMaskWidthX ){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2+iMaskWidthX && ix < iGapWidthX/2+iMaskWidthX+iTrenchWidthX ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2+iMaskWidthX+iTrenchWidthX && ix < iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX ){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX && ix < iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX*2 ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX*2 && ix < iGapWidthX/2+iMaskWidthX*3+iTrenchWidthX*2 ){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2+iMaskWidthX*3+iTrenchWidthX*2 && ix < iGapWidthX/2*3+iMaskWidthX*3+iTrenchWidthX*2 ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*3+iTrenchWidthX*2 && ix < iGapWidthX/2*3+iMaskWidthX*4+iTrenchWidthX*2){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*4+iTrenchWidthX*2 && ix < iGapWidthX/2*3+iMaskWidthX*4+iTrenchWidthX*3 ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*4+iTrenchWidthX*3 && ix < iGapWidthX/2*3+iMaskWidthX*5+iTrenchWidthX*3){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*5+iTrenchWidthX*3 && ix < iGapWidthX/2*3+iMaskWidthX*5+iTrenchWidthX*4 ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*5+iTrenchWidthX*4 && ix < iGapWidthX/2*3+iMaskWidthX*6+iTrenchWidthX*4){
                                                        C1.setStatus(itag, iMaskStat, iNumMask);
                                                }else if ( ix >= iGapWidthX/2*3+iMaskWidthX*6+iTrenchWidthX*4 && ix < iGapWidthX/2*4+iMaskWidthX*6+iTrenchWidthX*4){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }
                                        }
                                }
                        }
                }
        }else if ( meshfile == "single_finite_trench" ){
                C1.initial(Nx, Ny, Nz, Lx, Ly, Lz);
                for(int iz = 0; iz < Nz; iz++){
                        for(int iy = 0; iy < Ny; iy++){
                                for(int ix = 0; ix < Nx; ix++){
                                        int itag =  ix + ( iy + iz*Ny )*Nx;
                                        if(  iz < iSubstrateThickZ ){
                                                //--setting for substrate
                                                C1.setStatus(itag, iSubstrateStat, iNumMaterial);     //--Si has 8 atoms per unit cell with the dimension of 0.54 nm
                                        }else if (  iz >= iSubstrateThickZ && iz < (iSubstrateThickZ + iMaskThickZ)  ) {
                                                if( iy < iGapWidthY/2){
                                                        C1.setStatus(itag, iVacuumStat, 0);
                                                }else if ( iy >= iGapWidthY/2 && iy < iGapWidthY/2+iMaskLength ){
                                                        if(  ix < iGapWidthX/2 ){
                                                                //--setting for trench
                                                                C1.setStatus(itag, iVacuumStat, 0);
                                                        }else if ( ix >= iGapWidthX/2 && ix < iGapWidthX/2+iMaskWidthX ){
                                                                //--setting for mask
                                                                C1.setStatus(itag, iMaskStat, iNumMask);
                                                        }else if ( ix >= iGapWidthX/2+iMaskWidthX && ix < iGapWidthX/2+iMaskWidthX+iTrenchWidthX ){
                                                                //--setting for trench
                                                                C1.setStatus(itag, iVacuumStat, 0);
                                                        }else if ( ix >= iGapWidthX/2+iMaskWidthX+iTrenchWidthX && ix < iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX ){
                                                                //--setting for mask
                                                                C1.setStatus(itag, iMaskStat, iNumMask);
                                                        }else if ( ix >= iGapWidthX/2+iMaskWidthX*2+iTrenchWidthX && ix < iGapWidthX+iMaskWidthX*2+iTrenchWidthX){
                                                                //--setting for trench
                                                                C1.setStatus(itag, iVacuumStat, 0);
                                                        }
                                                }else if( iy >= iGapWidthY/2+iMaskLength && iy < iGapWidthY+iMaskLength ){
                                                        C1.setStatus(itag, iVacuumStat, 0);
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
                         }else{
                                 for( int i = 0; i < POINT_DATA; i++){
                                         C1.dNumSiClxs[i][0] = iNumMaterial;
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
                         }else{
                                 for( int i = 0; i < POINT_DATA; i++){
                                         C1.dNumSiClxs[i][1] = 0;
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
                         }else{
                                 for( int i = 0; i < POINT_DATA; i++){
                                         C1.dNumSiClxs[i][2] = 0;
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
                         }else{
                                 for( int i = 0; i < POINT_DATA; i++){
                                         C1.dNumSiClxs[i][3] = 0;
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
        double SpeedcutoffThermalParicle [ThermalParticleTypes]; //--Array to store pre-calculated cutoff speed
        SpeedcutoffThermalParicle[iClRadicalType] = CalculateSpeedCutoff(  MassChlorine, Temperature  );
        SpeedcutoffThermalParicle[iSigType]              = CalculateSpeedCutoff(  MassSilicon, Temperature  );
        SpeedcutoffThermalParicle[iSiClgType]          = CalculateSpeedCutoff(  MassSilicon+MassChlorine, Temperature  );
        SpeedcutoffThermalParicle[iSiCl2gType]       = CalculateSpeedCutoff(  MassSilicon+2*MassChlorine, Temperature  );
        SpeedcutoffThermalParicle[iSiCl3gType]       = CalculateSpeedCutoff(  MassSilicon+3*MassChlorine, Temperature  );


         /* Pre-calculation for particle generation probability */
        double TotalFlux = dClRadicalFlux+dClIonFlux+dCl2IonFlux+dArIonType;
        vector<double> GenerationProbIncidentParticle;
        GenerationProbIncidentParticle.push_back(dClRadicalFlux/TotalFlux);
        GenerationProbIncidentParticle.push_back(dClIonFlux/TotalFlux);
        GenerationProbIncidentParticle.push_back(dCl2IonFlux/TotalFlux);
        GenerationProbIncidentParticle.push_back(dArIonType/TotalFlux);




        /*Pre-calculation of particle size and initialization of six point particle*/
        double SixPointParticle [6][3] = {  {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1} };
        double ParticleSizeX, ParticleSizeY, ParticleSizeZ;
        ParticleSizeX = ParticleSizeFactor*dx;
        ParticleSizeY = ParticleSizeFactor*dy;
        ParticleSizeZ = ParticleSizeFactor*dz;


        /*Pre-calculation for surface sites*/
        int SurfaceSearchingIndex [SurfaceSearchingNumber][3];
        for (int k = 0; k < SurfaceSearchingRange; k++){
                for (int j = 0; j < SurfaceSearchingRange; j++){
                        for(int i = 0; i < SurfaceSearchingRange; i++){
                                SurfaceSearchingIndex[i+(j+k*SurfaceSearchingRange )*SurfaceSearchingRange ][0] = i - SurfaceSearchingRadius;
                                SurfaceSearchingIndex[i+(j+k*SurfaceSearchingRange )*SurfaceSearchingRange ][1] = j - SurfaceSearchingRadius;
                                SurfaceSearchingIndex[i+(j+k*SurfaceSearchingRange )*SurfaceSearchingRange ][2] = k -SurfaceSearchingRadius;
                        }
                }
        }


        /*Pre-calculation of total particle number, time interval,  and output frequency*/
        double FluxArea = Lx*Ly;
        double RealTimeInterval = 1/(TotalFlux*1E4)/FluxArea;  //--convert unit from cm^-2 to m^-2
        int TotalParticle;
        cout  << "The area of incidence (A) = " << FluxArea << " m^2" << endl;
        cout << "The total flux = "<< TotalFlux*1E4 << " m^-2 s^-1" << endl;
        //cout << "The number of atoms in a cell (Ns) = " << iNumMaterial << endl;
        cout << "The real time interval between steps (dt) = " << "1/(total flux * area of incidence) = " << RealTimeInterval << endl;
        TotalParticle = int(TotalRealTime/RealTimeInterval);
        cout << "The total number of steps in the simulation will be : " << TotalParticle << endl;
        int FileIndex = 0;
        int frequency = int(TotalParticle/OutputFileNumber);
        int particleNumber = 0;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                     PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+OutputFilename, append, FileIndex );




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
                P1.setInitialType( GenerationProbIncidentParticle );  //--choose a particular type of particle
                double rdd1, rdd2, rdd3;                                //--random number for selecting ion energy and angle
                int ReactionExecution = 0;                              //--determine if a reaction happen
                int EmittedParticle = 0;                                   //--index for emitted particle such as SiClx(g)
                int ReflectedParticle = 0;                              //--index for original reflected particle such as Ar*, Cl*
                int ReactionIndex = 0;                                 //--index for a particular reaction formula
                double normSurfaceNormal [3];                              //--normalized surface normal vector (x, y, z)
                double normReflectedVelocity [3];                            //--normalized reflected velocity vector (x, y, z)
                double GrazingAngle = 0;                               //--angle between surface and velocity
                double IncidentAngle = 0;                              //--angle between normal and velocity
                int old_iPos [3];                                       //--to record the iPos (x, y, z) before propogation
                double old_dPos [3];                                    //--to record the dPos (x, y, z) before propagation
                int itag;                                               //--itag of particle center
                int old_itag;
                int itag_six_point [6] ;                                //--itag of particle's six neighboring points
                double dPos_six_point [6][3] ;                          //--dPos (x, y, z) of particle's six neighboring points
                int iPos_six_point [6][3] ;                             //--iPos (x, y, z) of particle's siz neighboring points
                int CountPointInSolid;                               //--to count how many point of a seven-point molecule is on solid cell
                double normEmittedNormal [3];                         //--normalized velocity for emitted particle
                double dPos_EmittedParticle [3];                           //--dPos for emitted particle
                int iPos_EmittedParticle [3];                              //--iPos for emitted particle
                double prob_of_energy [10]= {0.0};  //--used in ion enhanced reaction, the energy-dependent reaction probability
                double prob_of_angle [10]= {0.0};  //--used in ion enhanced reaction, the angle-dependent reaction probability
                double accerlation [3];

                //--Count total number of particle and write to a file in a vtk format
                #pragma omp critical
                {
                        particleNumber++;

                        if(  particleNumber%frequency == 0  ){
                                FileIndex++;
                                if (  FileIndex == OutputFileNumber  ){
                                        PRINT_SI = true;
                                        PRINT_SICl = true;
                                        PRINT_SICl2 = true;
                                        PRINT_SICl3 = true;
                                }
                                write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                                                PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+OutputFilename, append, FileIndex );
                        }
                }


                if (P1.ParticleType == 0){
                        continue;
                }else if ( P1.ParticleType == iClRadicalType ){
                        P1.mass = MassChlorine;
                        P1.speed = P1.setInitialSpeed( Temperature,   MassChlorine,   SpeedcutoffThermalParicle[iClRadicalType]  );
                        P1.energy = 0.5 * P1.mass * P1.speed * P1.speed; //--unit: Joule
                        P1.theta = unif(generator)*PI/2+PI/2;  //--randomly generate a value between PI/2 to PI for the polar angle
                        if ( NEUTRAL_THETA_SCALING == true){

                                P1.theta = (P1.theta-PI)*NeutralThetaScalingFactor+PI;  //--unit: radian
                        }
                }else if ( P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType ){
                        rdd1 = unif(generator);
                        rdd2 = unif(generator);
                        rdd3 = unif(generator);
                        if ( P1.ParticleType == iClIonType){
                                P1.mass = MassChlorine;
                                for(int i = 0;    i < cumulativeflux_ClIon.size()-1;     i++){
                                        if  (  rdd1 >= cumulativeflux_ClIon[i] && rdd1 < cumulativeflux_ClIon[i+1]  ){  //--read cumulative flux data for Cl ion
                                                P1.theta =(    IonAngle[i] + rdd2*(  IonAngle[i+1] - IonAngle[i]  )    )*PI/180+PI;  //--unit: radian
                                                P1.energy = (    IonEnergy[i] + rdd3*(  IonEnergy[i+1] - IonEnergy[i]  )    )/Joule_to_eV;  //--unit: Joule
                                                break;
                                        }
                                }
                        }else if ( P1.ParticleType == iCl2IonType){
                                P1.mass = MassChlorine*2;
                                for(int i = 0;    i < cumulativeflux_Cl2Ion.size()-1;     i++){
                                        if  (  rdd1 >= cumulativeflux_Cl2Ion[i] && rdd1 < cumulativeflux_Cl2Ion[i+1]  ){  //--read cumulative flux data for Cl2 ion
                                                P1.theta =(    IonAngle[i] + rdd2*(  IonAngle[i+1] - IonAngle[i]  )    )*PI/180+PI;  //--unit: radian
                                                P1.energy = (    IonEnergy[i] + rdd3*(  IonEnergy[i+1] - IonEnergy[i]  )    )/Joule_to_eV;  //--unit: Joule
                                                break;
                                        }
                                }
                        }else if (  P1.ParticleType == iArIonType){
                                P1.mass = MassArgon;
                                for(int i = 0;    i < cumulativeflux_ArIon.size()-1;     i++){  //--read cumulative flux data for Ar ion
                                        if  (  rdd1 >= cumulativeflux_ArIon[i] && rdd1 < cumulativeflux_ArIon[i+1]  ){
                                                P1.theta =(    IonAngle[i] + rdd2*(  IonAngle[i+1] - IonAngle[i])    )*PI/180+PI;  //--unit: radian
                                                P1.energy = (    IonEnergy[i] + rdd3*(  IonEnergy[i+1] - IonEnergy[i]  )    )/Joule_to_eV;  //--unit: Joule
                                                break;
                                        }
                                }
                        }
                        if ( ION_THETA_SCALING == true){
                                P1.theta = (P1.theta-PI)*IonThetaScalingFactor+PI; //--unit: radian
                        }
                        //cout << "P1.theta of ion = " << P1.theta << endl;
                        //cin.get();
                        P1.speed = sqrt(2*P1.energy/P1.mass);
                }


                P1.phi = unif(generator)*2*PI; //--unit: radian
                P1.setInitialPosition(Lx, Ly, Lz, dx, dy, dz);
                P1.PropagationTimeInterval = dx/P1.speed;
                P1.Vel[X_dir] = P1.speed*sin(P1.theta)*cos(P1.phi);
                P1.Vel[Y_dir] = P1.speed*sin(P1.theta)*sin(P1.phi);
                P1.Vel[Z_dir] = P1.speed*cos(P1.theta);


                for(int indexTimeStep = 0; indexTimeStep < PropagationTimestepNumber; indexTimeStep++){

                        //cout << "indexTimeStep = " << indexTimeStep << endl;
                        //cout << "particle type = " << P1.ParticleType << endl;

                        //--check if particle has been deactivated
                        if(P1.ParticleType == 0){

                                if (EmittedParticle == 0){
                                        break;
                                }else{
                                        P1.ParticleType = EmittedParticle;
                                        if (P1.ParticleType == iSiClgType){
                                                P1.mass = MassSilicon+MassChlorine;
                                        }else if ( P1.ParticleType == iSiCl2gType){
                                                P1.mass = MassSilicon+2*MassChlorine;
                                        }else if ( P1.ParticleType == iSiCl3gType){
                                                P1.mass = MassSilicon+3*MassChlorine;
                                        }else if ( P1.ParticleType == iSiCl4gType){
                                                EmittedParticle = 0;
                                                break;
                                        }

                                        P1.dPos[X_dir] = dPos_EmittedParticle[X_dir];
                                        P1.dPos[Y_dir] = dPos_EmittedParticle[Y_dir];
                                        P1.dPos[Z_dir] = dPos_EmittedParticle[Z_dir];
                                        P1.iPos[X_dir] = iPos_EmittedParticle[X_dir];
                                        P1.iPos[Y_dir] = iPos_EmittedParticle[Y_dir];
                                        P1.iPos[Z_dir] = iPos_EmittedParticle[Z_dir];
                                        P1.speed = P1.setInitialSpeed( Temperature,   P1.mass,   SpeedcutoffThermalParicle[P1.ParticleType]  );
                                        P1.energy = 0.5*P1.mass*P1.speed*P1.speed;
                                        P1.theta = P1.setReemittedTheta(ReemissionCosineLawPower); //--theta is with respect to surface normal
                                        P1.phi = unif(generator)*2*PI;
                                        //--modify particle reflected velocity by new speed, theta, and phi
                                        P1.ReemittedWithNewDirection(normEmittedNormal, P1.speed, P1.theta, P1.phi) ;
                                        P1.PropagationTimeInterval = dx/P1.speed;
                                        EmittedParticle = 0;
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
                        if (P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType){
                                //cout << "particle type" << P1.ParticleType << endl;
                                accerlation[X_dir] = C1.ElectricForce[old_itag][X_dir] / P1.mass;
                                accerlation[Y_dir] = C1.ElectricForce[old_itag][Y_dir] / P1.mass;
                                accerlation[Z_dir] = C1.ElectricForce[old_itag][Z_dir] / P1.mass;
                                //cout << "accerlation[X_dir]  " << accerlation[X_dir]<< endl;
                                //cin.get();
                        }else{
                                accerlation[X_dir] = 0;
                                accerlation[Y_dir] = 0;
                                accerlation[Z_dir] = 0;
                        }

                        P1.Vel[X_dir]       -= accerlation[X_dir] * P1.PropagationTimeInterval;
                        P1.Vel[Y_dir]       -= accerlation[Y_dir] * P1.PropagationTimeInterval;
                        P1.Vel[Z_dir]       -= accerlation[Z_dir] * P1.PropagationTimeInterval;

                        P1.dPos[X_dir]   += P1.Vel[X_dir]*P1.PropagationTimeInterval-0.5*accerlation[X_dir]*pow(P1.PropagationTimeInterval, 2);
                        P1.dPos[Y_dir]   += P1.Vel[Y_dir]*P1.PropagationTimeInterval-0.5*accerlation[Y_dir]*pow(P1.PropagationTimeInterval, 2);
                        P1.dPos[Z_dir]   += P1.Vel[Z_dir]*P1.PropagationTimeInterval-0.5*accerlation[Z_dir]*pow(P1.PropagationTimeInterval, 2);

                        P1.iPos[X_dir]     = floor(P1.dPos[X_dir]/dx);
                        P1.iPos[Y_dir]     = floor(P1.dPos[Y_dir]/dy);
                        P1.iPos[Z_dir]     = floor(P1.dPos[Z_dir]/dz);
                       //cin.get();
                        //--boundary condition
                        if (BoundaryCondition == "reflective"){
                                if(  P1.iPos[X_dir] <  0 ){
                                        P1.dPos[X_dir] = P1.dPos[X_dir]*(-1);
                                        P1.iPos[X_dir] = floor(P1.dPos[X_dir]/dx);
                                        P1.Vel[X_dir] = P1.Vel[X_dir]*(-1);
                                }else if( P1.iPos[X_dir] >= Nx ){
                                        P1.dPos[X_dir] = 2*Lx - P1.dPos[X_dir];
                                        P1.iPos[X_dir] = floor(P1.dPos[X_dir]/dx);
                                        P1.Vel[X_dir] = P1.Vel[X_dir]*(-1);
                                }
                                if( P1.iPos[Y_dir] < 0 ){
                                        P1.dPos[Y_dir] = P1.dPos[Y_dir]*(-1);
                                        P1.iPos[Y_dir] = floor(P1.dPos[Y_dir]/dy);
                                        P1.Vel[Y_dir] = P1.Vel[Y_dir]*(-1);
                                }else if( P1.iPos[Y_dir] >= Ny ){
                                        P1.dPos[Y_dir] = 2*Ly - P1.dPos[Y_dir];
                                        P1.iPos[Y_dir] = floor(P1.dPos[Y_dir]/dy);
                                        P1.Vel[Y_dir] = P1.Vel[Y_dir]*(-1);
                                }
                        }else if (BoundaryCondition == "periodic"){
                                if(P1.iPos[X_dir] < 0){
                                        P1.dPos[X_dir] += Lx;
                                        P1.iPos[X_dir] += Nx;
                                }else if(P1.iPos[X_dir] >= Nx){
                                        P1.dPos[X_dir] -= Lx;
                                        P1.iPos[X_dir] -= Nx;
                                }
                                if(P1.iPos[Y_dir] < 0){
                                        P1.dPos[Y_dir] += Ly;
                                        P1.iPos[Y_dir] += Ny;
                                }else if(P1.iPos[Y_dir] >= Ny){
                                        P1.dPos[Y_dir] -= Ly;
                                        P1.iPos[Y_dir] -= Ny;
                                }
                        }else{
                                if (  P1.iPos[X_dir] <  0 ||  P1.iPos[X_dir] >= Nx ){
                                        P1.ParticleType = 0;
                                        continue;
                                }
                                if (   P1.iPos[Y_dir] < 0 ||  P1.iPos[Y_dir] >= Ny ){
                                        P1.ParticleType = 0;
                                        continue;
                                }
                        }

                        if(P1.iPos[Z_dir] < 0 || P1.iPos[Z_dir] >= Nz){
                                P1.ParticleType = 0;
                                continue;
                        }else{
                                itag =  P1.iPos[X_dir] + ( P1.iPos[Y_dir] + P1.iPos[Z_dir]*Ny )*Nx;
                        }


                        //--Calculate particle neighboring six points (particle size edge)
                        for ( int i = 0; i < 6 ; i++){
                                dPos_six_point[i][X_dir]   = P1.dPos[X_dir] + SixPointParticle[i][X_dir]*ParticleSizeX;
                                dPos_six_point[i][Y_dir]   = P1.dPos[Y_dir] + SixPointParticle[i][Y_dir]*ParticleSizeY;
                                dPos_six_point[i][Z_dir]   = P1.dPos[Z_dir] + SixPointParticle[i][Z_dir]*ParticleSizeZ;
                                iPos_six_point[i][X_dir]     = floor(dPos_six_point[i][X_dir]/dx  );
                                iPos_six_point[i][Y_dir]     = floor( dPos_six_point[i][Y_dir]/dy  );
                                iPos_six_point[i][Z_dir]     = floor( dPos_six_point[i][Z_dir]/dz );


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
                                        itag_six_point[i] =  iPos_six_point[i][X_dir] + ( iPos_six_point[i][Y_dir] + iPos_six_point[i][Z_dir]*Ny )*Nx;
                                }
                        }


                        //--count how many point for a seven-point molecule is on solid
                        CountPointInSolid =0;
                        for (int i = 0; i < 6 ; i++){
                                if (C1.iStatus[itag_six_point[i]] == iSubstrateStat || C1.iStatus[itag_six_point[i]] == iMaskStat){
                                        CountPointInSolid++;
                                }
                        }
                        if (C1.iStatus[itag] == iSubstrateStat || C1.iStatus[itag] == iMaskStat ){
                                CountPointInSolid++;
                        }



                        //--check if a particle collide on solid cell
                        if (CountPointInSolid > 1){
                                P1.dPos[X_dir] = old_dPos[X_dir];
                                P1.dPos[Y_dir] = old_dPos[Y_dir];
                                P1.dPos[Z_dir] = old_dPos[Z_dir];
                                P1.PropagationTimeInterval = 0.5 * P1.PropagationTimeInterval;
                        }else if ( CountPointInSolid == 1){
                                for (int i = 0; i < 6 ; i++){
                                        if (C1.iStatus[itag_six_point[i]] == iSubstrateStat || C1.iStatus[itag_six_point[i]] == iMaskStat){
                                                itag = itag_six_point[i];
                                        }
                                }
                                C1.surface_normal(SurfaceSearchingIndex, SurfaceSearchingNumber , itag, P1.iPos, P1.Vel, normSurfaceNormal,
                                                                           normReflectedVelocity, &GrazingAngle, &IncidentAngle );
                                /*
                                double n_dot_vi;
                                C1.CalculateSurfaceNormal(SurfaceSearchingIndex, SurfaceSearchingNumber, itag, P1.iPos, normSurfaceNormal );
                                C1.CalculateReflectedVelocity(P1.Vel, normSurfaceNormal, normReflectedVelocity, &n_dot_vi);

                                if ( acos(n_dot_vi)*180/PI  >  90 ){
                                        IncidentAngle = 180 - acos(n_dot_vi)*180/PI;
                                }else{
                                        IncidentAngle = acos(n_dot_vi)*180/PI;
                                }
                                GrazingAngle =  90 - IncidentAngle;
                                */
                                /*
                                if (P1.iPos[Z_dir] < 950 && P1.ParticleType == iClIonType){
                                        cout << "particle position = " << P1.iPos[X_dir] << " " << P1.iPos[Y_dir] << " " << P1.iPos[Z_dir] << endl;
                                        cout << "normSurfaceNormal = " << normSurfaceNormal[X_dir] << " " << normSurfaceNormal[Y_dir] << " " << normSurfaceNormal[Z_dir]<<endl;
                                        cout << "normIncidentV = " << P1.Vel[X_dir]/P1.speed << " " << P1.Vel[Y_dir]/P1.speed << " " << P1.Vel[Z_dir]/P1.speed << endl;
                                        cout << "normReflectedV = " << normReflectedVelocity[X_dir] << " " << normReflectedVelocity[Y_dir] << " " << normReflectedVelocity[Z_dir]<<endl;
                                        cout << "IncidentAngle = " << IncidentAngle << endl;
                                        cout << "GrazinAngel = " << GrazingAngle << endl;
                                        cin>>n_dot_vi;
                                }
                                */


                                if (C1.iStatus[itag] == iMaskStat){

                                        if (P1.ParticleType == iClRadicalType){
                                                P1.Vel[X_dir] = P1.speed*normReflectedVelocity[X_dir];
                                                P1.Vel[Y_dir] = P1.speed*normReflectedVelocity[Y_dir];
                                                P1.Vel[Z_dir] = P1.speed*normReflectedVelocity[Z_dir];
                                        }


                                        if (P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType ){
                                                break; //--for test
                                                C1.IonMaskReaction(PhysSputterProb, itag, P1.energy*Joule_to_eV, IncidentAngle, &ReactionExecution  );
                                                //--modify particle energy, speed, and reflected velocity
                                                P1.ReflectedWithNewEnergy(normReflectedVelocity, &GrazingAngle, &epsilon_0, &epsilon_s, &theta_0, &gamma_0);
                                        }


                                        if (P1.ParticleType == iSiClgType || P1.ParticleType == iSiCl2gType || P1.ParticleType == iSiCl3gType ){
                                                break; //--for test
                                        }


                                        if ( P1.speed == 0){
                                                P1.ParticleType = 0;
                                        }else{
                                                P1.PropagationTimeInterval = dx/P1.speed;
                                        }
                                        continue;


                                }else if ( C1.iStatus[itag] == iSubstrateStat ){

                                        if( P1.ParticleType == iClRadicalType ){
                                                C1.ClRadicalReaction(p0_ClRadicalReaction, itag, iNumMaterial, &ReactionExecution, &ReactionIndex);
                                                if (ReactionExecution == 0){
                                                        P1.speed = P1.setInitialSpeed( Temperature,   P1.mass,   SpeedcutoffThermalParicle[P1.ParticleType]  );
                                                        P1.energy = 0.5*P1.mass*P1.speed*P1.speed;
                                                        P1.theta = P1.setReemittedTheta(ReemissionCosineLawPower); //--theta is with respect to surface normal
                                                        P1.phi = unif(generator)*2*PI;
                                                        //--modify particle reflected velocity by new speed, theta, and phi
                                                        P1.ReemittedWithNewDirection(normSurfaceNormal, P1.speed, P1.theta, P1.phi) ;
                                                        P1.PropagationTimeInterval = dx/P1.speed;
                                                        continue;
                                                }else if (ReactionExecution == 1){
                                                        P1.ParticleType = 0;
                                                        continue;
                                                }

                                        }else if ( P1.ParticleType == iClIonType || P1.ParticleType == iCl2IonType || P1.ParticleType == iArIonType){

                                                if (P1.ParticleType == iClIonType){

                                                        C1.ClIonReaction(Eth_ClIonReaction, &E0_ClIonReaction, p0_ClIonReaction, type_ClIonReaction,
                                                                                               PhysSputterProb, ChemSputterProb, itag, iNumMaterial, P1.energy*Joule_to_eV, IncidentAngle,
                                                                                               &ReactionExecution, &ReflectedParticle, &EmittedParticle, &ReactionIndex);
                                                        /*
                                                        calc_prob_of_energy(Eth_ClIonReaction , &E0_ClIonReaction, P1.energy*Joule_to_eV, prob_of_energy);
                                                        calc_prob_of_angle( type_ClIonReaction, PhysSputterProb, ChemSputterProb, &IncidentAngle, prob_of_angle);
                                                        C1.ClIonReaction(itag, p0_ClIonReaction, prob_of_energy, prob_of_angle, &ReactionExecution, &ReflectedParticle,
                                                                                              &EmittedParticle, &ReactionIndex);
                                                        */
                                                }else if (P1.ParticleType == iCl2IonType){

                                                        C1.Cl2IonReaction(Eth_Cl2IonReaction, &E0_Cl2IonReaction, p0_Cl2IonReaction, type_Cl2IonReaction,
                                                                                                  PhysSputterProb, ChemSputterProb, itag, iNumMaterial, P1.energy*Joule_to_eV, IncidentAngle,
                                                                                                  &ReactionExecution, &ReflectedParticle, &EmittedParticle, &ReactionIndex);
                                                        /*
                                                        calc_prob_of_energy(Eth_Cl2IonReaction,  &E0_Cl2IonReaction, P1.energy*Joule_to_eV, prob_of_energy);
                                                        calc_prob_of_angle(type_Cl2IonReaction, PhysSputterProb, ChemSputterProb, &IncidentAngle, prob_of_angle);
                                                        C1.Cl2IonReaction(itag, p0_Cl2IonReaction, prob_of_energy, prob_of_angle, &ReactionExecution, &ReflectedParticle,
                                                                                                &EmittedParticle, &ReactionIndex);
                                                        */
                                                }else if (P1.ParticleType == iArIonType){
                                                        C1.ArIonReaction(Eth_ArIonReaction, &E0_ArIonReaction, p0_ArIonReaction, type_ArIonReaction,
                                                                                               PhysSputterProb, ChemSputterProb, itag, iNumMaterial, P1.energy*Joule_to_eV, IncidentAngle,
                                                                                               &ReactionExecution, &ReflectedParticle, &EmittedParticle, &ReactionIndex);
                                                        /*
                                                        calc_prob_of_energy(Eth_ArIonReaction,  &E0_ArIonReaction, P1.energy*Joule_to_eV, prob_of_energy);
                                                        calc_prob_of_angle( type_ArIonReaction, PhysSputterProb, ChemSputterProb, &IncidentAngle, prob_of_angle);
                                                        C1.ArIonReaction(itag, p0_ArIonReaction, prob_of_energy, prob_of_angle, &ReactionExecution, &ReflectedParticle,
                                                                                              &EmittedParticle, &ReactionIndex);
                                                        */
                                                }

                                                if (ReactionExecution == 1 && EmittedParticle != 0){
                                                        dPos_EmittedParticle[X_dir] = P1.dPos[X_dir];
                                                        dPos_EmittedParticle[Y_dir] = P1.dPos[Y_dir];
                                                        dPos_EmittedParticle[Z_dir] = P1.dPos[Z_dir];
                                                        iPos_EmittedParticle[X_dir] = P1.iPos[X_dir];
                                                        iPos_EmittedParticle[Y_dir] = P1.iPos[Y_dir];
                                                        iPos_EmittedParticle[Z_dir] = P1.iPos[Z_dir];
                                                        normEmittedNormal[X_dir] = normSurfaceNormal[X_dir];
                                                        normEmittedNormal[Y_dir] = normSurfaceNormal[Y_dir];
                                                        normEmittedNormal[Z_dir] = normSurfaceNormal[Z_dir];
                                                }

                                                P1.ParticleType = ReflectedParticle;
                                                if (P1.ParticleType == 0){
                                                        continue;
                                                }else{
                                                        if (P1.ParticleType == iClIonType){
                                                                P1.mass = MassChlorine;
                                                        }else if (P1.ParticleType == iCl2IonType){
                                                                P1.mass = MassChlorine*2;
                                                        }else if (P1.ParticleType == iArIonType){
                                                                P1.mass = MassArgon;
                                                        }
                                                        //--modify particle energy, speed, and reflected velocity
                                                        P1.ReflectedWithNewEnergy(normReflectedVelocity, &GrazingAngle, &epsilon_0, &epsilon_s, &theta_0, &gamma_0);
                                                        if ( P1.speed == 0){
                                                                P1.ParticleType = 0;
                                                        }else{
                                                                P1.PropagationTimeInterval = dx/P1.speed;
                                                        }
                                                        continue;
                                                }

                                        }else if ( P1.ParticleType == iSiClgType || P1.ParticleType == iSiCl2gType || P1.ParticleType == iSiCl3gType){
                                                C1.redeposition(p0_redeposition, P1.ParticleType, itag, &ReactionExecution, &ReactionIndex);

                                                if (ReactionExecution == 0){
                                                        P1.speed = P1.setInitialSpeed( Temperature, P1.mass, SpeedcutoffThermalParicle[P1.ParticleType]  );
                                                        P1.energy = 0.5*P1.mass*P1.speed*P1.speed;
                                                        P1.theta = P1.setReemittedTheta(ReemissionCosineLawPower); //--theta is with respect to surface normal
                                                        P1.phi = unif(generator)*2*PI;
                                                        //--modify particle reflected velocity by new speed, theta, and phi
                                                        P1.ReemittedWithNewDirection(normSurfaceNormal, P1.speed, P1.theta, P1.phi) ;
                                                        P1.PropagationTimeInterval = dx/P1.speed;
                                                        continue;
                                                }else if (ReactionExecution == 1){
                                                        P1.ParticleType = 0;
                                                        continue;
                                                }

                                        }//--End of emittion particle ( if clause)
                                }//--End of substrate (if clause)
                        }//--End of count_point_on_solid (if clause)
                }//--End of particle propagation loop (for loop)
        }//--End of particle generation loop (for loop)
        cout << "Total number of files generated : " << FileIndex <<endl;
        cout << "Total Particle Number : " << particleNumber << endl;
        FileIndex++;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                     PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+OutputFilename, append, FileIndex );

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
