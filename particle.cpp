#include "cell.h"
#include "etching.h"
#include "rand.h"
#include "random"

using namespace std;


void cell::Generation(int TotalParticleInAFile, double Temperature, double* CumuProb,
                                    vector<double>& cumulativeflux_ClIon, vector<double>& cumulativeflux_Cl2Ion, vector<double>& cumulativeflux_ArIon,
                                    vector<double>& IonEnergy, vector<double>& IonAngle,
                                    bool NEUTRAL_THETA_SCALING, bool ION_THETA_SCALING,
                                    double NeutralThetaScalingFactor, double IonThetaScalingFactor,
                                    int* ParticleType, double* mass, double** dPos, double** Vel, double* energy )
{

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0, 1.0);

        for(int iParticle = 0; iParticle < TotalParticleInAFile; iParticle++){

                double randParticle = ran3(&iRandTag);
                double randPosX = ran3(&iRandTag);
                double randPosY = ran3(&iRandTag);
                double randVel1 = ran3(&iRandTag);
                double randVel2 = ran3(&iRandTag);
                double randVel3 = ran3(&iRandTag);
                double theta, phi, speed;

                if (randParticle >= 0 && randParticle < CumuProb[0]){

                        ParticleType[iParticle] = iClRadicalType;
                        mass[iParticle] = MassChlorine;
                        dPos[iParticle][X_dir] = randPosX*cell::dDimLength[X_dir];
                        dPos[iParticle][Y_dir] = randPosY*cell::dDimLength[Y_dir];
                        dPos[iParticle][Z_dir] = cell::dDimLength[Z_dir]-cell::dDimLength_per_cell[Z_dir]*0.5;
                        do{

                                double v_mp = sqrt(2*BoltzmannConstant_SI*Temperature/mass[iParticle]);  //--the most probable velocity = sqrt(2kT/m)
                                Vel[iParticle][X_dir] = v_mp*distribution(generator);
                                Vel[iParticle][Y_dir] = v_mp*distribution(generator);
                                Vel[iParticle][Z_dir] = v_mp*distribution(generator);
                                speed = sqrt(Vel[iParticle][X_dir]*Vel[iParticle][X_dir]+Vel[iParticle][Y_dir]*Vel[iParticle][Y_dir]+Vel[iParticle][Z_dir]*Vel[iParticle][Z_dir]  );
                        }while( speed == 0);
                        if (Vel[iParticle][Z_dir] > 0)    Vel[iParticle][Z_dir] *= (-1);
                        if ( NEUTRAL_THETA_SCALING == true  && theta > 0){
                                theta = (theta-PI)*NeutralThetaScalingFactor+PI;  //--unit: radian
                                phi = randVel3*2*PI;
                                Vel[iParticle][X_dir] = speed*sin(theta)*cos(phi);
                                Vel[iParticle][Y_dir] = speed*sin(theta)*sin(phi);
                                Vel[iParticle][Z_dir] = speed*cos(theta);
                        }
                        energy[iParticle] = 0.5*mass[iParticle]*speed*speed; //--unit: Joule
                }else if ( randParticle >= CumuProb[0] && randParticle < CumuProb[1]){

                        ParticleType[iParticle] = iClIonType;
                        mass[iParticle] = MassChlorine;
                        dPos[iParticle][X_dir] = randPosX*cell::dDimLength[X_dir];
                        dPos[iParticle][Y_dir] = randPosY*cell::dDimLength[Y_dir];
                        dPos[iParticle][Z_dir] = cell::dDimLength[Z_dir]-cell::dDimLength_per_cell[Z_dir]*0.5;
                        for(int grid = 0;    grid < cumulativeflux_ClIon.size()-1;     grid++){
                                if  (  randVel1 >= cumulativeflux_ClIon[grid] && randVel1 < cumulativeflux_ClIon[grid+1]  ){  //--read cumulative flux data for Cl ion
                                        theta =(    IonAngle[grid] + randVel2*(  IonAngle[grid+1] - IonAngle[grid]  )    )*PI/180+PI;  //--unit: radian
                                        energy[iParticle] = (    IonEnergy[grid] + randVel2*(  IonEnergy[grid+1] - IonEnergy[grid]  )    )/Joule_to_eV;  //--unit: Joule
                                        break;
                                }
                        }
                        if ( ION_THETA_SCALING == true && theta > 0){
                                theta = (theta-PI)*IonThetaScalingFactor+PI; //--unit: radian
                        }
                        speed = sqrt(2*energy[iParticle]/mass[iParticle]);
                        phi = randVel3*2*PI; //--unit: radian
                        Vel[iParticle][X_dir] = speed*sin(theta)*cos(phi);
                        Vel[iParticle][Y_dir] = speed*sin(theta)*sin(phi);
                        Vel[iParticle][Z_dir] = speed*cos(theta);

                }else if ( randParticle >= CumuProb[1] && randParticle < CumuProb[2]){

                        ParticleType[iParticle] = iCl2IonType;
                        mass[iParticle] = MassChlorine*2;
                        dPos[iParticle][X_dir] = randPosX*cell::dDimLength[X_dir];
                        dPos[iParticle][Y_dir] = randPosY*cell::dDimLength[Y_dir];
                        dPos[iParticle][Z_dir] = cell::dDimLength[Z_dir]-cell::dDimLength_per_cell[Z_dir]*0.5;
                        for(int grid = 0;    grid < cumulativeflux_Cl2Ion.size()-1;     grid++){
                                if  (  randVel1 >= cumulativeflux_Cl2Ion[grid] && randVel1 < cumulativeflux_Cl2Ion[grid+1]  ){  //--read cumulative flux data for Cl2 ion
                                        theta = (    IonAngle[grid] + randVel2*(  IonAngle[grid+1] - IonAngle[grid]  )    )*PI/180+PI;  //--unit: radian
                                        energy[iParticle] = (    IonEnergy[grid] + randVel2*(  IonEnergy[grid+1] - IonEnergy[grid]  )    )/Joule_to_eV;  //--unit: Joule
                                        break;
                                }
                        }
                        if ( ION_THETA_SCALING == true && theta > 0){
                                theta = (theta - PI)*IonThetaScalingFactor+PI; //--unit: radian
                        }
                        speed = sqrt(2*energy[iParticle]/mass[iParticle]);
                        phi = randVel3*2*PI; //--unit: radian
                        Vel[iParticle][X_dir] = speed*sin(theta)*cos(phi);
                        Vel[iParticle][Y_dir] = speed*sin(theta)*sin(phi);
                        Vel[iParticle][Z_dir] = speed*cos(theta);
                }else if ( randParticle >= CumuProb[2] && randParticle < CumuProb[3]){

                        ParticleType[iParticle] = iArIonType;
                        mass[iParticle] = MassArgon;
                        dPos[iParticle][X_dir] = randPosX*cell::dDimLength[X_dir];
                        dPos[iParticle][Y_dir] = randPosY*cell::dDimLength[Y_dir];
                        dPos[iParticle][Z_dir] = cell::dDimLength[Z_dir]-cell::dDimLength_per_cell[Z_dir]*0.5;
                        for(int grid = 0;    grid < cumulativeflux_ArIon.size()-1;     grid++){  //--read cumulative flux data for Ar ion
                                if  (  randVel1 >= cumulativeflux_ArIon[grid] && randVel1 < cumulativeflux_ArIon[grid+1]  ){
                                        theta =(    IonAngle[grid] + randVel2*(  IonAngle[grid+1] - IonAngle[grid])    )*PI/180+PI;  //--unit: radian
                                        energy[iParticle] = (    IonEnergy[grid] + randVel2*(  IonEnergy[grid+1] - IonEnergy[grid]  )    )/Joule_to_eV;  //--unit: Joule
                                        break;
                                }
                        }
                        if ( ION_THETA_SCALING == true){
                                theta = (theta-PI)*IonThetaScalingFactor+PI; //--unit: radian
                        }
                        speed = sqrt(2*energy[iParticle]/mass[iParticle]);
                        phi = randVel3*2*PI; //--unit: radian
                        Vel[iParticle][X_dir] = speed*sin(theta)*cos(phi);
                        Vel[iParticle][Y_dir] = speed*sin(theta)*sin(phi);
                        Vel[iParticle][Z_dir] = speed*cos(theta);

                }else{
                        ParticleType[iParticle] = 0;
                }
        }
}



void cell::Propagation(int TotalParticleInAFile, double Temperature, int PropagationTimestepNumber, double ReemissionCosineLawPower,
                                      double ParticleSizeFactor, int SurfaceSearchingRadius,
                                      int* ParticleType, double* mass, double** dPos, double** Vel, double* energy )
{


        /*Pre-calculation of particle size and initialization of six point particle*/
        double indexParticleEdge [6][3] = {  {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1} };
        double ParticleSizeX, ParticleSizeY, ParticleSizeZ;
        ParticleSizeX = ParticleSizeFactor*cell::dDimLength_per_cell[X_dir];
        ParticleSizeY = ParticleSizeFactor*cell::dDimLength_per_cell[Y_dir];
        ParticleSizeZ = ParticleSizeFactor*cell::dDimLength_per_cell[Z_dir];

        /*Pre-calculation for surface sites*/
        int SurfaceSearchingRange = 2*SurfaceSearchingRadius+1;
        int SurfaceSearchingNumber = pow(  SurfaceSearchingRange, 3);
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

        for (int iParticle = 0;   iParticle < TotalParticleInAFile;   iParticle++){

                cout << "iParticle = " << iParticle << endl;
                //--To record the itag and position of the particle center
                int iTag;                                                                                  //--itag of particle center
                int iLastTag;                                                                         //--itag of particle cneter of the previous step
                double dLastPos [3];                                                       //--dPos (x, y, z) of the previous step
                int iLastPos [3];                                                                  //--iPos (x, y, z) of the previous step

                //--To record the itag and position of the particle edge
                int iParticleEdgeTag [6] ;                                               //--itag of particle's neighboring six edge points
                double dParticleEdgePos [6][3] ;                              //--dPos (x, y, z) of particle's neighboring six edge points
                int iParticleEdgePos [6][3] ;                                         //--iPos (x, y, z) of particle's six neighboring edge points
                int CountPointInSolid = 0;                                            //--to count how many point of a seven-point molecule is on solid cell

                //--To record the normalized normal vector and reflected vector
                double normSurfaceNormal [3] = {0.0};               //--normalized surface normal vector (x, y, z)
                double normReflectedVelocity [3] = {0.0};          //--normalized reflected velocity vector (x, y, z)
                double GrazingAngle = 0.0;                                         //--angle between surface and velocity
                double IncidentAngle = 0.0;                                         //--angle between normal and velocity

                //-Calculate spped, PropagationTimeInterval and iPos for each particle
                double speed = sqrt( pow(Vel[iParticle][X_dir], 2.0)+pow(Vel[iParticle][Y_dir], 2.0)+pow(Vel[iParticle][Z_dir], 2.0) );
                double PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/speed;
                int iPos [3];
                iPos[X_dir] = int(dPos[iParticle][X_dir]/cell::dDimLength_per_cell[X_dir]);
                iPos[Y_dir] = int(dPos[iParticle][Y_dir]/cell::dDimLength_per_cell[Y_dir]);
                iPos[Z_dir] = int(dPos[iParticle][Z_dir]/cell::dDimLength_per_cell[Z_dir]);



                for(int indexTimeStep = 0; indexTimeStep < PropagationTimestepNumber; indexTimeStep++){

                        //--Check if particle has been deactivated
                        if(ParticleType[iParticle] == 0)    break;

                        //--Set iLastPos and dLastPos to be current iPos and dPos
                        iLastPos[X_dir] = iPos[X_dir];
                        iLastPos[Y_dir] = iPos[Y_dir];
                        iLastPos[Z_dir] = iPos[Z_dir];
                        dLastPos[X_dir] =dPos[iParticle][X_dir] ;
                        dLastPos[Y_dir] =dPos[iParticle][Y_dir] ;
                        dLastPos[Z_dir] =dPos[iParticle][Z_dir] ;
                        iLastTag = iLastPos[X_dir] + ( iLastPos[Y_dir] + iLastPos[Z_dir]*cell::iDimSize[Y_dir] )*cell::iDimSize[X_dir];

                        //--Propagation
                        dPos[iParticle][X_dir]   += Vel[iParticle][X_dir]*PropagationTimeInterval ;
                        dPos[iParticle][Y_dir]   += Vel[iParticle][Y_dir]*PropagationTimeInterval;
                        dPos[iParticle][Z_dir]   += Vel[iParticle][Z_dir]*PropagationTimeInterval;
                        iPos[X_dir]     = floor(  dPos[iParticle][X_dir]/cell::dDimLength_per_cell[X_dir]);
                        iPos[Y_dir]     = floor(  dPos[iParticle][Y_dir]/cell::dDimLength_per_cell[Y_dir]);
                        iPos[Z_dir]     = floor(  dPos[iParticle][Z_dir]/cell::dDimLength_per_cell[Z_dir]);

                        cell::BoundaryMapping(iPos, dPos[iParticle], &ParticleType[iParticle]);
                        cell::CalcParticleEdge(dPos[iParticle], &ParticleType[iParticle], iParticleEdgeTag, indexParticleEdge,
                                                             ParticleSizeX, ParticleSizeY, ParticleSizeZ, dParticleEdgePos, iParticleEdgePos  ) ;
                        iTag =  iPos[X_dir] + ( iPos[Y_dir] + iPos[Z_dir]*iDimSize[Y_dir] )*iDimSize[X_dir];
                        if (ParticleType[iParticle] == 0)   continue;


                        //--Count how many point for a seven-point molecule is on solid
                        CountPointInSolid=0;
                        if (cell::iStatus[iParticleEdgeTag[0]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[0]] == iMaskStat)    CountPointInSolid++;
                        if (cell::iStatus[iParticleEdgeTag[1]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[1]] == iMaskStat)     CountPointInSolid++;
                        if (cell::iStatus[iParticleEdgeTag[2]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[2]] == iMaskStat)    CountPointInSolid++;
                        if (cell::iStatus[iParticleEdgeTag[3]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[3]] == iMaskStat)    CountPointInSolid++;
                        if (cell::iStatus[iParticleEdgeTag[4]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[4]] == iMaskStat)    CountPointInSolid++;
                        if (cell::iStatus[iParticleEdgeTag[5]] == iSubstrateStat || cell::iStatus[iParticleEdgeTag[5]] == iMaskStat)    CountPointInSolid++;
                        if (cell::iStatus[iTag] == iSubstrateStat                                   || cell::iStatus[iTag] == iMaskStat )                                    CountPointInSolid++;


                        //--check if a particle collide on solid cell
                        if (CountPointInSolid  > 1){
                                dPos[iParticle][X_dir] = dLastPos[X_dir];
                                dPos[iParticle][Y_dir] = dLastPos[Y_dir];
                                dPos[iParticle][Z_dir] = dLastPos[Z_dir];
                                PropagationTimeInterval = 0.5 * PropagationTimeInterval;
                        }else if ( CountPointInSolid == 1){
                                for (int i = 0; i < 6 ; i++){
                                        if ( iStatus[iParticleEdgeTag[i]] == iSubstrateStat || iStatus[iParticleEdgeTag[i]] == iMaskStat){
                                                iLastTag = iTag;
                                                iTag = iParticleEdgeTag[i];
                                                break;
                                        }
                                }
                                CalcSurfaceNormal(SurfaceSearchingIndex, SurfaceSearchingNumber, iTag, iPos, Vel[iParticle],
                                normSurfaceNormal, normReflectedVelocity, &GrazingAngle, &IncidentAngle );
                                SurfaceReaction( iTag, iLastTag, GrazingAngle, IncidentAngle, &ParticleType[iParticle], &mass[iParticle], energy[iParticle]);

                                if(ParticleType[iParticle] < 6 ){//--iClRadical, iSiClg, iSiCl2g, iSiCl3g
                                        if ( iStatus[iTag] == iMaskStat){
                                                ElasticReflection(normReflectedVelocity, speed, Vel[iParticle] );
                                                PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/speed;
                                                continue;
                                        }else if ( iStatus[iTag] == iSubstrateStat){
                                                DiffusiveReflection(Temperature, normSurfaceNormal, ReemissionCosineLawPower, &energy[iParticle], mass[iParticle], Vel[iParticle] );
                                                PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/speed;
                                                continue;
                                        }
                                }else if ( ParticleType[iParticle] > 6){//--iClIon, iCl2Ion, iArIon
                                        if ( iStatus[iTag] == iMaskStat){
                                                ElasticReflection(normReflectedVelocity, speed, Vel[iParticle] );
                                                PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/speed;
                                                continue;
                                        }else if ( iStatus[iTag] == iSubstrateStat){
                                                InelasticReflection(normReflectedVelocity, &GrazingAngle, &energy[iParticle], mass[iParticle], Vel[iParticle] );
                                                if ( speed == 0){
                                                        ParticleType[iParticle] = 0;
                                                }else{
                                                        PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/speed;
                                                }
                                                continue;
                                        }
                                }else if  ( ParticleType[iParticle] == 0){
                                        continue;
                                }
                        }
                }
        }
}



void cell::BoundaryMapping(int* iPos, double* dPos, int* ParticleType){

        //--periodic boundary condition
        if(  iPos[X_dir] < 0  ){
                dPos[X_dir] += cell::dDimLength[X_dir];
                iPos[X_dir] += cell::iDimSize[X_dir];
        }else if(  iPos[X_dir] >= cell::iDimSize[X_dir]  ){
                dPos[X_dir] -= cell::dDimLength[X_dir];
                iPos[X_dir] -= cell::iDimSize[X_dir];
        }
        if(  iPos[Y_dir] < 0  ){
                dPos[Y_dir] += cell::dDimLength[Y_dir];
                iPos[Y_dir] += cell::iDimSize[Y_dir];
        }else if(  iPos[Y_dir] >= cell::iDimSize[Y_dir]){
                dPos[Y_dir] -= cell::dDimLength[Y_dir];
                iPos[Y_dir] -= cell::iDimSize[Y_dir];
        }
        if(  iPos[Z_dir] < 0 || iPos[Z_dir] >= cell::iDimSize[Z_dir]){
                *ParticleType = 0;
        }
}



void cell::CalcParticleEdge(double* dPos, int* ParticleType, int* iParticleEdgeTag,  double indexParticleEdge [][3],
                                                double ParticleSizeX, double ParticleSizeY, double ParticleSizeZ, double dParticleEdgePos [][3], int iParticleEdgePos [][3]  )
{

        //--Calculate particle neighboring six points (particle size edge)
        for ( int i = 0; i < 6 ; i++){
                dParticleEdgePos[i][X_dir]   = dPos[X_dir] + indexParticleEdge[i][X_dir]*ParticleSizeX;
                dParticleEdgePos[i][Y_dir]   = dPos[Y_dir] + indexParticleEdge[i][Y_dir]*ParticleSizeY;
                dParticleEdgePos[i][Z_dir]   = dPos[Z_dir] + indexParticleEdge[i][Z_dir]*ParticleSizeZ;
                iParticleEdgePos[i][X_dir]     = floor(dParticleEdgePos[i][X_dir]/cell::dDimLength_per_cell[X_dir]  );
                iParticleEdgePos[i][Y_dir]     = floor(dParticleEdgePos[i][Y_dir]/cell::dDimLength_per_cell[Y_dir]  );
                iParticleEdgePos[i][Z_dir]     = floor(dParticleEdgePos[i][Z_dir]/cell::dDimLength_per_cell[Z_dir] );

                if(iParticleEdgePos[i][X_dir] < 0){
                        dParticleEdgePos[i][X_dir] = dParticleEdgePos[i][X_dir] + cell::dDimLength[X_dir];
                        iParticleEdgePos[i][X_dir] = iParticleEdgePos[i][X_dir] + cell::iDimSize[X_dir];
                }else if(  iParticleEdgePos[i][X_dir] >= cell::iDimSize[X_dir]){
                        dParticleEdgePos[i][X_dir] = dParticleEdgePos[i][X_dir] - cell::dDimLength[X_dir];
                        iParticleEdgePos[i][X_dir] = iParticleEdgePos[i][X_dir] - cell::iDimSize[X_dir];
                }
                if( iParticleEdgePos[i][Y_dir] < 0){
                        dParticleEdgePos[i][Y_dir] = dParticleEdgePos[i][Y_dir] + cell::dDimLength[Y_dir];
                        iParticleEdgePos[i][Y_dir] = iParticleEdgePos[i][Y_dir] + cell::iDimSize[Y_dir];
                }else if( iParticleEdgePos[i][Y_dir] >= cell::iDimSize[Y_dir]){
                        dParticleEdgePos[i][Y_dir] = dParticleEdgePos[i][Y_dir] - cell::dDimLength[Y_dir];
                        iParticleEdgePos[i][Y_dir] = iParticleEdgePos[i][Y_dir] - cell::iDimSize[Y_dir];
                }
                if(iParticleEdgePos[i][Z_dir] < 0 || iParticleEdgePos[i][Z_dir] >= cell::iDimSize[Z_dir]){
                        *ParticleType = 0;
                        continue;
                }else{
                        iParticleEdgeTag[i] =  iParticleEdgePos[i][X_dir] +
                                                                      ( iParticleEdgePos[i][Y_dir] + iParticleEdgePos[i][Z_dir]*cell::iDimSize[Y_dir] )*cell::iDimSize[X_dir];
                }
        }
}




void cell::InelasticReflection(double* normReflectedVelocity, double* GrazingAngle, double* energy , double mass, double* Vel)
{

        double EnergyDependentScalingFactor;
        double AngleDependentScalingFactor;
        double Energy = *energy*Joule_to_eV;
        double epsilon_0 = 0.0;  //--unit: eV
        double epsilon_s = 50;  //--unit: eV
        double theta_0 = 30;  //--unit: degree
        double gamma_0 = 0.85;  //--a scaling factor, no unit
        double speed;
        if ( Energy < epsilon_0 ){
                EnergyDependentScalingFactor = 0;
        }else if ( Energy >= epsilon_0 && Energy <= epsilon_s ){
                EnergyDependentScalingFactor = (Energy - epsilon_0)/(epsilon_s - epsilon_0);
        }else if (  Energy > epsilon_s){
                EnergyDependentScalingFactor = 1;
        }

        if( *GrazingAngle  < theta_0){
                AngleDependentScalingFactor = (theta_0 - *GrazingAngle)/(theta_0);
        }else if (*GrazingAngle >= theta_0){
                AngleDependentScalingFactor = 0.0;
        }

        *energy = *energy * gamma_0 * AngleDependentScalingFactor * EnergyDependentScalingFactor; //--updated energy
        speed = sqrt(2* *energy/mass);                                                                                                                                              //--updated speed
        Vel[X_dir] = speed*normReflectedVelocity[X_dir];                                                                                                     //--updated x component of velocity
        Vel[Y_dir] = speed*normReflectedVelocity[Y_dir];                                                                                                     //--updated y component of velocity
        Vel[Z_dir] = speed*normReflectedVelocity[Z_dir];                                                                                                     //--updated z component of
}



void cell::ElasticReflection(double* normReflectedVelocity, double speed, double* Vel)
{
        Vel[X_dir] = speed*normReflectedVelocity[X_dir];
        Vel[Y_dir] = speed*normReflectedVelocity[Y_dir];
        Vel[Z_dir] = speed*normReflectedVelocity[Z_dir];

}


void cell::DiffusiveReflection(double Temperature, double* normSurfaceNormal, double ReemissionCosineLawPower,
                                                    double* energy, double mass, double* Vel)
{

        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0, 1.0);
        double randTheta = ran3(&iRandTag);
        double randPhi = ran3(&iRandTag);
        double speed;
        double theta = acos(   pow( randTheta, 1/(ReemissionCosineLawPower+1) )   );
        double phi = randPhi*2*PI;
        double alpha;  //--azimusal angle of surface normal vector
        double beta;   //--polar angle of surface normal vector
        double sqrt_of_Nxsquare_plus_Nysquare = sqrt(  pow(normSurfaceNormal[X_dir], 2.0)+pow(normSurfaceNormal[Y_dir], 2.0)  );
        if(  sqrt_of_Nxsquare_plus_Nysquare == 0 ){
                alpha = 0;
                beta = 0;
        }else{
                alpha = acos(normSurfaceNormal[X_dir] / sqrt_of_Nxsquare_plus_Nysquare );
                beta = acos(normSurfaceNormal[Z_dir]);
        }
        do{
                double v_mp = sqrt(2*BoltzmannConstant_SI*Temperature/mass);  //--the most probable velocity = sqrt(2kT/m)
                Vel[X_dir] = v_mp*distribution(generator);
                Vel[Y_dir] = v_mp*distribution(generator);
                Vel[Z_dir] = v_mp*distribution(generator);
                speed = sqrt(Vel[X_dir]*Vel[X_dir]+Vel[Y_dir]*Vel[Y_dir]+Vel[Z_dir]*Vel[Z_dir]);
        }while(speed == 0);
        *energy = 0.5*mass*speed*speed;


        //--modify particle reflected velocity by new speed, theta, and phi
        //--Use the formula that I derived
        Vel[X_dir] = speed*(   cos(alpha)*cos(beta)*sin(theta)*cos(phi)+cos(alpha)*sin(beta)*cos(theta)-sin(alpha)*sin(theta)*sin(phi)   );
        Vel[Y_dir] = speed*(   sin(alpha)*cos(beta)*sin(theta)*cos(phi)+sin(alpha)*sin(beta)*cos(theta)+cos(alpha)*sin(theta)*sin(phi)   );
        Vel[Z_dir] = speed*( -sin(beta)*sin(theta)*cos(phi)+cos(beta)*cos(theta)  );
}






