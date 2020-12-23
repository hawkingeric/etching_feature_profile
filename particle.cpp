#include "particle.h"
#include "etching.h"
#include "rand.h"
#include "random"

using namespace std;



//--input: normReflectedVelocity, GrazingAngle, epsilon_0, epsilon_s, theta_0, gamma_0, E, speed, velocity;  output: ReflectedVelocity
#pragma acc routine seq
void ReflectedWithNewEnergy(double* normReflectedVelocity, double* GrazingAngle, int indexParticle,
                                                                   double* energy, double* speed, double* mass, double* Vel){
        double EnergyDependentScalingFactor;
        double AngleDependentScalingFactor;
        double Energy = energy[indexParticle]*Joule_to_eV;
        int X_dir = 0, Y_dir = 1, Z_dir = 2;
        double epsilon_0 = 0.0;  //--unit: eV
        double epsilon_s = 50;  //--unit: eV
        double theta_0 = 30;  //--unit: degree
        double gamma_0 = 0.85;  //--a scaling factor, no unit

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

        energy[indexParticle] *= gamma_0*AngleDependentScalingFactor*EnergyDependentScalingFactor; //--updated energy
        *speed = sqrt(2*energy[indexParticle]/mass[indexParticle]); //--updated speed
        Vel[X_dir] = *speed*normReflectedVelocity[X_dir];  //--updated x component of velocity
        Vel[Y_dir] = *speed*normReflectedVelocity[Y_dir];  //--updated y component of velocity
        Vel[Z_dir] = *speed*normReflectedVelocity[Z_dir];  //--updated z component of velocity
}

#pragma acc routine seq
void ReemittedWithNewDirection(double* normSurfaceNormal, double speed, double theta, double phi, double* Vel){
        double alpha;  //--azimusal angle of surface normal vector
        double beta;   //--polar angle of surface normal vector
        int X_dir = 0, Y_dir = 1, Z_dir = 2;
        double sqrt_of_Nxsquare_plus_Nysquare = sqrt(  pow(normSurfaceNormal[X_dir], 2.0)+pow(normSurfaceNormal[Y_dir], 2.0)  );
        if(  sqrt_of_Nxsquare_plus_Nysquare == 0 ){
                alpha = 0;
                beta = 0;
        }else{
                alpha = acos(normSurfaceNormal[X_dir] / sqrt_of_Nxsquare_plus_Nysquare );
                beta = acos(normSurfaceNormal[Z_dir]);
        }

        //--Use the formula that I derived
        Vel[X_dir] = speed*(   cos(alpha)*cos(beta)*sin(theta)*cos(phi)+cos(alpha)*sin(beta)*cos(theta)-sin(alpha)*sin(theta)*sin(phi)   );
        Vel[Y_dir] = speed*(   sin(alpha)*cos(beta)*sin(theta)*cos(phi)+sin(alpha)*sin(beta)*cos(theta)+cos(alpha)*sin(theta)*sin(phi)   );
        Vel[Z_dir] = speed*( -sin(beta)*sin(theta)*cos(phi)+cos(beta)*cos(theta)  );
}


void Generation(int TotalParticleInAFile, double Temperature, double dx, double dy, double dz, double Lx, double Ly, double Lz, double* CumuProb,
                                    vector<double>& cumulativeflux_ClIon, vector<double>& cumulativeflux_Cl2Ion, vector<double>& cumulativeflux_ArIon,
                                    vector<double>& IonEnergy, vector<double>& IonAngle,
                                    bool NEUTRAL_THETA_SCALING, bool ION_THETA_SCALING,
                                    double NeutralThetaScalingFactor, double IonThetaScalingFactor,
                                    int* ParticleType, double* mass, double** Vel, double* energy, double** dPos)
{
        int X_dir = 0, Y_dir = 1, Z_dir = 2;
        int iClRadicalType = 1, iSigType = 2, iSiClgType = 3, iSiCl2gType = 4, iSiCl3gType = 5;
        int iClIonType = 7, iCl2IonType = 8, iArIonType = 9;
        double MassChlorine = 35.45*1.66053904E-27;  // unit: kg
        double MassArgon = 39.948*1.66053904E-27;  // unit: kg
        double MassSilicon = 28.0855*1.66053904E-27;  // unit: kg
        double MassElectron = 9.10938356E-31; // unit: kg

        for(int iParticle = 0; iParticle < TotalParticleInAFile; iParticle++){
                double randParticle = ran3(&iRandTag);
                double randPosX = ran3(&iRandTag);
                double randPosY = ran3(&iRandTag);
                double randVel1 = ran3(&iRandTag);
                double randVel2 = ran3(&iRandTag);
                double randVel3 = ran3(&iRandTag);
                double theta, phi, speed;
                cout << randParticle << " " << randPosX << " " << randPosY << " " << randVel1 << " " << randVel2 << " " << randVel3 << endl;

                if (randParticle >= 0 && randParticle < CumuProb[0]){

                        ParticleType[iParticle] = iClRadicalType;
                        mass[iParticle] = MassChlorine;
                        dPos[iParticle][X_dir] = randPosX*Lx;
                        dPos[iParticle][Y_dir] = randPosY*Ly;
                        dPos[iParticle][Z_dir] = Lz-dz*0.5;
                        while( speed == 0 ){
                                std::default_random_engine generator;
                                std::normal_distribution<double> distribution(0.0, 1.0);
                                double v_mp = sqrt(2*BoltzmannConstant_SI*Temperature/mass[iParticle]);  //--the most probable velocity = sqrt(2kT/m)
                                Vel[iParticle][X_dir] = v_mp*distribution(generator);
                                Vel[iParticle][Y_dir] = v_mp*distribution(generator);
                                Vel[iParticle][Z_dir] = v_mp*distribution(generator);
                        }
                        if (Vel[iParticle][Z_dir] > 0)    Vel[iParticle][Z_dir]*(-1);
                        speed = sqrt(Vel[iParticle][X_dir]*Vel[iParticle][X_dir]+Vel[iParticle][Y_dir]*Vel[iParticle][Y_dir]+Vel[iParticle][Z_dir]*Vel[iParticle][Z_dir]  );
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
                        dPos[iParticle][X_dir] = randPosX*Lx;
                        dPos[iParticle][Y_dir] = randPosY*Ly;
                        dPos[iParticle][Z_dir] = Lz-dz*0.5;
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
                        dPos[iParticle][X_dir] = randPosX*Lx;
                        dPos[iParticle][Y_dir] = randPosY*Ly;
                        dPos[iParticle][Z_dir] = Lz-dz*0.5;
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
                        dPos[iParticle][X_dir] = randPosX*Lx;
                        dPos[iParticle][Y_dir] = randPosY*Ly;
                        dPos[iParticle][Z_dir] = Lz-dz*0.5;
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



void Propagation(int TotalParticleInAFile, double Temperature, int PropagationTimestepNumber, double ReemissionCosineLawPower,
                                      int* ParticleType, double* mass, double** dPos, double* energy, double** Vel,
                                      double dx, double dy, double dz, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz,
                                      double SixPointParticle [][3], double ParticleSizeX, double ParticleSizeY, double ParticleSizeZ,
                                      int* iStatus, int** iID_NBR, double* dNumMaterial, double** dNumSiClxs, double** dNumSiClxg,
                                      int SurfaceSearchingIndex [][3] , int SurfaceSearchingNumber)
{


        int X_dir = 0, Y_dir = 1, Z_dir = 2;
        int iClRadicalType = 1, iSigType = 2, iSiClgType = 3, iSiCl2gType = 4, iSiCl3gType = 5;
        int iClIonType = 7, iCl2IonType = 8, iArIonType = 9;
        int iVacuumStat = 0, iSubstrateStat = 1, iMaskStat = 2;
        double MassChlorine = 35.45*1.66053904E-27;  // unit: kg
        double MassArgon = 39.948*1.66053904E-27;  // unit: kg
        double MassSilicon = 28.0855*1.66053904E-27;  // unit: kg
        double MassElectron = 9.10938356E-31; // unit: kg

        for (int indexParticle = 0;   indexParticle < TotalParticleInAFile;   indexParticle++){


                double normSurfaceNormal [3] = {0.0};                              //--normalized surface normal vector (x, y, z)
                double normReflectedVelocity [3] = {0.0};                            //--normalized reflected velocity vector (x, y, z)
                double GrazingAngle = 0.0;                               //--angle between surface and velocity
                double IncidentAngle = 0.0;                              //--angle between normal and velocity
                int old_iPos [3];                                       //--to record the iPos (x, y, z) before propogation
                double old_dPos [3];                                    //--to record the dPos (x, y, z) before propagation
                int itag;                                               //--itag of particle center
                int old_itag;
                int itag_six_point [6] ;                                //--itag of particle's six neighboring points
                double dPos_six_point [6][3] ;                          //--dPos (x, y, z) of particle's six neighboring points
                int iPos_six_point [6][3] ;                             //--iPos (x, y, z) of particle's siz neighboring points
                int CountPointInSolid;                               //--to count how many point of a seven-point molecule is on solid cell
                double speed = sqrt( pow(Vel[indexParticle][X_dir], 2.0)+pow(Vel[indexParticle][Y_dir], 2.0)+pow(Vel[indexParticle][Z_dir], 2.0) );
                double PropagationTimeInterval = dx/speed;
                double randTheta = ran3(&iRandTag);
                double randPhi = ran3(&iRandTag);
                double theta, phi;
                int iPos [3];
                iPos[X_dir] = int(dPos[indexParticle][X_dir]/dx);
                iPos[Y_dir] = int(dPos[indexParticle][Y_dir]/dy);
                iPos[Z_dir] = int(dPos[indexParticle][Z_dir]/dz);



                for(int indexTimeStep = 0; indexTimeStep < PropagationTimestepNumber; indexTimeStep++){

                        //--check if particle has been deactivated
                        if(ParticleType[indexParticle] == 0){
                                break;
                        }

                        //--Give old_iPos and old_dPos values
                        old_iPos[X_dir] = iPos[X_dir];
                        old_iPos[Y_dir] = iPos[Y_dir];
                        old_iPos[Z_dir] = iPos[Z_dir];
                        old_dPos[X_dir] =dPos[indexParticle][X_dir] ;
                        old_dPos[Y_dir] =dPos[indexParticle][Y_dir] ;
                        old_dPos[Z_dir] =dPos[indexParticle][Z_dir] ;
                        old_itag = old_iPos[X_dir] + ( old_iPos[Y_dir] + old_iPos[Z_dir]*Ny )*Nx;


                        //--propagation
                        dPos[indexParticle][X_dir]   += Vel[indexParticle][X_dir]*PropagationTimeInterval ;
                        dPos[indexParticle][Y_dir]   += Vel[indexParticle][Y_dir]*PropagationTimeInterval;
                        dPos[indexParticle][Z_dir]   += Vel[indexParticle][Z_dir]*PropagationTimeInterval;
                        iPos[X_dir]     = floor(  dPos[indexParticle][X_dir]/dx);
                        iPos[Y_dir]     = floor(  dPos[indexParticle][Y_dir]/dy);
                        iPos[Z_dir]     = floor(  dPos[indexParticle][Z_dir]/dz);

                        //--boundary condition
                        if(  iPos[X_dir] < 0  ){
                                dPos[indexParticle][X_dir] += Lx;
                                iPos[X_dir] += Nx;
                        }else if(  iPos[X_dir] >= Nx  ){
                                dPos[indexParticle][X_dir] -= Lx;
                                iPos[X_dir] -= Nx;
                        }
                        if(  iPos[Y_dir] < 0  ){
                                dPos[indexParticle][Y_dir] += Ly;
                                iPos[Y_dir] += Ny;
                        }else if(  iPos[Y_dir] >= Ny){
                                dPos[indexParticle][Y_dir] -= Ly;
                                iPos[Y_dir] -= Ny;
                        }

                        if(  iPos[Z_dir] < 0 || iPos[Z_dir] >= Nz){
                                ParticleType[indexParticle] = 0;
                                continue;
                        }else{
                                itag =  iPos[X_dir] + ( iPos[Y_dir] + iPos[Z_dir]*Ny )*Nx;
                        }


                        //--Calculate particle neighboring six points (particle size edge)
                        for ( int i = 0; i < 6 ; i++){
                                dPos_six_point[i][X_dir]   = dPos[indexParticle][X_dir] + SixPointParticle[i][X_dir]*ParticleSizeX;
                                dPos_six_point[i][Y_dir]   = dPos[indexParticle][Y_dir] + SixPointParticle[i][Y_dir]*ParticleSizeY;
                                dPos_six_point[i][Z_dir]   = dPos[indexParticle][Z_dir] + SixPointParticle[i][Z_dir]*ParticleSizeZ;
                                iPos_six_point[i][X_dir]     = floor(dPos_six_point[i][X_dir]/dx  );
                                iPos_six_point[i][Y_dir]     = floor( dPos_six_point[i][Y_dir]/dy  );
                                iPos_six_point[i][Z_dir]     = floor( dPos_six_point[i][Z_dir]/dz );


                                if(iPos_six_point[i][X_dir] < 0){
                                        dPos_six_point[i][X_dir] = dPos_six_point[i][X_dir] + Lx;
                                        iPos_six_point[i][X_dir] = iPos_six_point[i][X_dir] + Nx;
                                }else if(  iPos_six_point[i][X_dir] >= Nx){
                                        dPos_six_point[i][X_dir] = dPos_six_point[i][X_dir] - Lx;
                                        iPos_six_point[i][X_dir] = iPos_six_point[i][X_dir] - Nx;
                                }
                                if( iPos_six_point[i][Y_dir] < 0){
                                        dPos_six_point[i][Y_dir] = dPos_six_point[i][Y_dir] + Ly;
                                        iPos_six_point[i][Y_dir] = iPos_six_point[i][Y_dir] + Ny;
                                }else if( iPos_six_point[i][Y_dir] >= Ny){
                                        dPos_six_point[i][Y_dir] = dPos_six_point[i][Y_dir] - Ly;
                                        iPos_six_point[i][Y_dir] = iPos_six_point[i][Y_dir] - Ny;
                                }
                                if(iPos_six_point[i][Z_dir] < 0 || iPos_six_point[i][Z_dir] >= Nz){
                                        ParticleType[indexParticle] = 0;
                                        continue;
                                }else{
                                        itag_six_point[i] =  iPos_six_point[i][X_dir] + ( iPos_six_point[i][Y_dir] + iPos_six_point[i][Z_dir]*Ny )*Nx;
                                }
                        }

                        //--count how many point for a seven-point molecule is on solid
                        CountPointInSolid=0;
                        for (  int i = 0; i < 6 ; i++){
                                if (iStatus[itag_six_point[i]] == iSubstrateStat || iStatus[itag_six_point[i]] == iMaskStat){
                                        CountPointInSolid++;
                                }
                        }
                        if (  iStatus[itag] == iSubstrateStat || iStatus[itag] == iMaskStat ){
                                CountPointInSolid++;
                        }

                        //--check if a particle collide on solid cell
                        if (CountPointInSolid  > 1){
                                dPos[indexParticle][X_dir] = old_dPos[X_dir];
                                dPos[indexParticle][Y_dir] = old_dPos[Y_dir];
                                dPos[indexParticle][Z_dir] = old_dPos[Z_dir];
                                PropagationTimeInterval = 0.5 * PropagationTimeInterval;
                        }else if ( CountPointInSolid == 1){
                                for (int i = 0; i < 6 ; i++){
                                        if ( iStatus[itag_six_point[i]] == iSubstrateStat || iStatus[itag_six_point[i]] == iMaskStat){
                                                old_itag = itag;
                                                itag = itag_six_point[i];
                                        }
                                }
                                surface_normal(SurfaceSearchingIndex, SurfaceSearchingNumber, itag, iPos, Vel[indexParticle],
                                                                    Nx, Ny, Nz, iStatus, iID_NBR, iVacuumStat, iSubstrateStat, iMaskStat,
                                                                    normSurfaceNormal, normReflectedVelocity, &GrazingAngle, &IncidentAngle );

                                surface_reaction( itag, old_itag, GrazingAngle, IncidentAngle, &ParticleType[indexParticle], &mass[indexParticle],
                                                                      &energy[indexParticle], iStatus, dNumMaterial, dNumSiClxs);

                                if(ParticleType[indexParticle] < 6 ){//--iClRadical, iSiClg, iSiCl2g, iSiCl3g
                                        if ( iStatus[itag] == iMaskStat){
                                                Vel[indexParticle][X_dir] = speed*normReflectedVelocity[X_dir];
                                                Vel[indexParticle][Y_dir] = speed*normReflectedVelocity[Y_dir];
                                                Vel[indexParticle][Z_dir] = speed*normReflectedVelocity[Z_dir];
                                                continue;
                                        }else{
                                                while( speed == 0 ){
                                                        std::default_random_engine generator;
                                                        std::normal_distribution<double> distribution(0.0, 1.0);
                                                        double v_mp = sqrt(2*BoltzmannConstant_SI*Temperature/mass[indexParticle]);  //--the most probable velocity = sqrt(2kT/m)
                                                        Vel[indexParticle][X_dir] = v_mp*distribution(generator);
                                                        Vel[indexParticle][Y_dir] = v_mp*distribution(generator);
                                                        Vel[indexParticle][Z_dir] = v_mp*distribution(generator);
                                                        speed = sqrt(Vel[indexParticle][X_dir]*Vel[indexParticle][X_dir]
                                                                                +Vel[indexParticle][Y_dir]*Vel[indexParticle][Y_dir]+Vel[indexParticle][Z_dir]*Vel[indexParticle][Z_dir]);
                                                }
                                                energy[indexParticle] = 0.5*mass[indexParticle]*speed*speed;
                                                theta = acos(   pow( randTheta, 1/(ReemissionCosineLawPower+1) )   );
                                                phi = randPhi*2*PI;
                                                //--modify particle reflected velocity by new speed, theta, and phi
                                                ReemittedWithNewDirection(normSurfaceNormal, speed, theta, phi, Vel[indexParticle] ) ;
                                                PropagationTimeInterval = dx/speed;
                                                continue;
                                        }
                                }else if ( ParticleType[indexParticle] > 6){//--iClIon, iCl2Ion, iArIon
                                        //--modify particle energy, speed, and reflected velocity
                                        ReflectedWithNewEnergy(normReflectedVelocity, &GrazingAngle, indexParticle, energy, &speed, mass, Vel[indexParticle]);
                                        if ( speed == 0){
                                                ParticleType[indexParticle] = 0;
                                        }else{
                                                PropagationTimeInterval = dx/speed;
                                        }
                                        continue;
                                }else if  ( ParticleType[indexParticle] == 0){
                                        continue;
                                }
                        }
                }//--End of particle propagation loop (for loop)
        }//--End of particle generation loop
}







