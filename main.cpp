#include <stdlib.h>
#include <iomanip>
#include <string>
#include "etching.h"
#include "cell.h"
#include "particle.h"
#include "json.hpp"
#include "rand.h"
#include "queue"
#include "openacc.h"
#include "math.h"
#include "cmath"
#define N 100
#define numBlocks 1
#define threadsPerBlock 100







using namespace std;
using json = nlohmann::json;
int X_dir = 0, Y_dir = 1, Z_dir = 2;
int iClRadicalType = 1, iSigType = 2, iSiClgType = 3, iSiCl2gType = 4, iSiCl3gType = 5;
int iClIonType = 7, iCl2IonType = 8, iArIonType = 9;
int iVacuumStat = 0, iSubstrateStat = 1, iMaskStat = 2;
double MassChlorine = 35.45*1.66053904E-27;  // unit: kg
double MassArgon = 39.948*1.66053904E-27;  // unit: kg
double MassSilicon = 28.0855*1.66053904E-27;  // unit: kg
double MassElectron = 9.10938356E-31; // unit: kg
int P_sputtering = 0, C_sputtering = 1;

int main(int argc, char* argv[])
{


        /* Declaration*/
        //--parameter read from mesh
        string inputJsonFile, tmpString, geometry, meshfile;
        double dx, dy, dz;
        int Nx, Ny, Nz, iSubstrateThickZ, iMaskThickZ, iTrenchWidthX, iMaskWidthX, iMaskLength, iGapWidthX, iGapWidthY,
               iNumMaterial, iNumMask;
        int SurfaceSearchingRadius;

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
                        ParticleSizeFactor, ReemissionCosineLawPower, TotalRealTime, FileTimeInterval;
        int PropagationTimestepNumber;

        //--parameters read from surface ractions
        vector<double> PhysSputterProb, ChemSputterProb;
        vector<double> p0_ClRadicalReaction, p0_redeposition, p0_ClIonReaction, p0_Cl2IonReaction, p0_ArIonReaction;
        vector<double> Eth_ClIonReaction, Eth_Cl2IonReaction, Eth_ArIonReaction;  //--Eth is the activation energy for ion-enhanced reactions
        double E0_ClIonReaction, E0_Cl2IonReaction, E0_ArIonReaction;  //--E0 is the reference energy
        vector<int> type_ClIonReaction, type_Cl2IonReaction, type_ArIonReaction ;  //--Type of reactions: either physical or chemical sputtering



        double Lx, Ly, Lz ;





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
                        FileTimeInterval = config["file_time_interval"];
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
                        TotalRealTime = config[ "total_real_time"];
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
                cout << "  cell_size                    =    " << dx << " meter" << endl ;
                cout << "  domain size x                =    " << Nx << endl ;
                cout << "  domain size y                =    " << Ny << endl ;
                cout << "  domain size z                =    " << Nz << endl ;
                cout << "  substrate_thickness          =    " << iSubstrateThickZ << endl;
                cout << "  mask_thickness               =    " << iMaskThickZ << endl;
                cout << "  trench_width                 =    " << iTrenchWidthX << endl;
                cout << "  mask_width                   =    " << iMaskWidthX << endl;
                cout << "  mask_length                  =    " << iMaskLength << endl;
                cout << "  gap_width_x                  =    " << iGapWidthX << endl;
                cout << "  gap_width_y                  =    " << iGapWidthY << endl;
                cout << "  atom_number_in_material      =    " << iNumMaterial << endl;
                cout << "  atom_number_in_mask          =    " << iNumMask << endl;
                cout << "  surface_searching_radius     =    " << SurfaceSearchingRadius << " number of cells" << endl;
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
                cout << "  file_time_interval           =    " << FileTimeInterval << " seconds" << endl;
                cout << endl;
                cout << "Simulation conditions:" << endl ;
                cout << "  total_real_time              =    "  << TotalRealTime << " secsonds" << endl;
                cout << "  propagation_timestep_number  =    " << PropagationTimestepNumber << endl ;
                cout << "  Boundary condition           =    " << BoundaryCondition << endl;
                cout << "  Temperature                  =    "  << Temperature << " Kelvin" << endl;
                cout << "  Cl radical flux              =    " << dClRadicalFlux << " cm^-2 s^-1" << endl;
	            cout << "  Cl ion flux                  =    " <<        dClIonFlux << " cm^-2 s^-1" << endl;
                cout << "  Cl2+ ion flux                =    " << dCl2IonFlux << " cm^-2 s^-1" << endl;
                cout <<  "  Ar ion flux                  =    " <<    dArIonType << " cm^-2 s^-1" <<endl;
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
                        //cout << data << endl;
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

                        //cout << data << endl;
                        if (data == "SCALARS dNumMaterial float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                array_NumMaterial.push_back(stoi(data));
                                        }
                                }
                         }

                         if (data == "SCALARS iStatus float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                status = stoi(data);
                                                if ( status == iVacuumStat || status == iSubstrateStat ){
                                                        C1.setStatus(i, status, array_NumMaterial[i]);
                                                }else if ( status == iMaskStat){
                                                        C1.setStatus(i, status, iNumMask);
                                                }
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSis float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][0] = stod(data)*iNumMaterial;
                                        }
                                }
                         }


                         if (data == "SCALARS dNumSiCls float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][1] = stod(data)*iNumMaterial;
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSiCl2s float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][2] = stod(data)*iNumMaterial;
                                        }
                                }
                         }

                         if (data == "SCALARS dNumSiCl3s float 1"  ){
                                getline(in, data);
                                if ( data == "LOOKUP_TABLE default" ){
                                        for( int i = 0; i < POINT_DATA; i++){
                                                getline(in, data);
                                                C1.dNumSiClxs[i][3] = stod(data)*iNumMaterial;
                                        }
                                }
                         }
	            }
	            in.close();

        }


        /* Pre-calculation of speed cutoffe */
        int ThermalParticleTypes = 6; //--Total number of thermal species, here we include Cl radical, Si(g), SiCl(g), SiCl2(g), and SiCl3(g)


         /* Pre-calculation for particle generation probability */
        double TotalFlux = dClRadicalFlux+dClIonFlux+dCl2IonFlux+dArIonType;
        double CumuProbIncidentParticle [4];
        CumuProbIncidentParticle[0] = dClRadicalFlux/TotalFlux;
        CumuProbIncidentParticle[1] = dClIonFlux/TotalFlux + CumuProbIncidentParticle[0];
        CumuProbIncidentParticle[2] = dCl2IonFlux/TotalFlux + CumuProbIncidentParticle[1];
        CumuProbIncidentParticle[3] = dArIonType/TotalFlux + CumuProbIncidentParticle[2];


        /*Pre-calculation of total particle number, time interval,  and output frequency*/
        double FluxArea = Lx*Ly;
        double material_density = 5.0*1E28; //--density of Si in unit of m^-3
        double ParticleTimeInterval = material_density*dx*dy*dz/(TotalFlux*1E4)/FluxArea;  //--convert flux unit from cm^-2 to m^-2
        int TotalParticleInAFile;
        cout  << "The area of incidence (A) = " << FluxArea << " m^2" << endl;
        cout << "The total flux (J) = "<< TotalFlux*1E4 << " m^-2 s^-1" << endl;
        cout << "The number of atoms in a pseudo-particle (Ns) = " << int(material_density*dx*dy*dz) << endl;
        cout << "The real time interval between each pseudo-particle particle (dt) = " << "Ns/(J*A) = " << ParticleTimeInterval << endl;
        int TotalOutputFile;
        TotalOutputFile = TotalRealTime/FileTimeInterval;
        TotalParticleInAFile = int(FileTimeInterval/ParticleTimeInterval);
        cout << "The total number of particles in a file_interval = " << TotalParticleInAFile << endl;
        int FileIndex = 0;
        //int frequency = int(TotalParticle/OutputFileNumber);
        int particleNumber = 0;
        write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                     PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+OutputFilename, append, FileIndex );


        ran3_ini(0);
        for ( int indexOutputFile = 1; indexOutputFile <= TotalOutputFile; indexOutputFile++){

                int* ParticleType = new int [TotalParticleInAFile];
                double* mass = new double [TotalParticleInAFile];
                double** dPos = new double* [TotalParticleInAFile];
                double** Vel = new double* [TotalParticleInAFile];
                double* energy = new double [TotalParticleInAFile];
                for(int i =0; i<TotalParticleInAFile; i++){
                        dPos[i] = new double [3];
                        Vel[i] = new double [3];
                }

                C1.Generation(TotalParticleInAFile, Temperature, CumuProbIncidentParticle,cumulativeflux_ClIon, cumulativeflux_Cl2Ion,
                                                cumulativeflux_ArIon, IonEnergy, IonAngle,NEUTRAL_THETA_SCALING, ION_THETA_SCALING,
                                                NeutralThetaScalingFactor, IonThetaScalingFactor, ParticleType, mass, dPos, Vel, energy);

                C1.Propagation(TotalParticleInAFile, Temperature, PropagationTimestepNumber, ReemissionCosineLawPower, ParticleSizeFactor,
                                                   SurfaceSearchingRadius, ParticleType, mass, dPos, Vel, energy);



                if (   indexOutputFile == TotalOutputFile  ){
                        PRINT_SI = true;
                        PRINT_SICl = true;
                        PRINT_SICl2 = true;
                        PRINT_SICl3 = true;
                }
                write_to_vtk(vtkDataType, Nx, Ny, Nz, dx, C1.iStatus, C1.dNumMaterial, C1.dNumMask, C1.dNumSiClxs, iNumMaterial,
                                    PRINT_SI, PRINT_SICl, PRINT_SICl2, PRINT_SICl3, directory+OutputFilename, append, indexOutputFile );
        }
        cout << "Total number of files generated : " << TotalOutputFile <<endl;
        cout << "Total Particle Number : " << particleNumber << endl;
}//--End of Program
