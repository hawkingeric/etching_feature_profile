#include "cell.h"
#include "etching.h"
#include "jacobi_eigenvalue.h"
#include "rand.h"
#include "random"

using namespace std;


void cell::Generation(double* Temperature, double* CumuProb, double* v_cut,
                                    vector<double> cumulativeflux_ClIon, vector<double> cumulativeflux_Cl2Ion, vector<double> cumulativeflux_ArIon,
                                    vector<double> IonEnergy, vector<double> IonAngle,
                                    bool NEUTRAL_THETA_SCALING, bool ION_THETA_SCALING,
                                    double NeutralThetaScalingFactor, double IonThetaScalingFactor,
                                    int* ParticleType, double* mass, double* dPos, int* iPos, double* Vel, double* speed, double* energy )
{
        double theta, phi;
        cell::GenerateParticleTypeMass(CumuProb, ParticleType, mass);
        cell::GeneratePosition(dPos,  iPos);
        if(*ParticleType == iClRadicalType){
                double ThetaScalingFactor = 1;
                cell::GenerateRadicalSpeed(Temperature, v_cut, mass, speed);
                cell::GenerateTheta(&theta);
                theta = (theta-PI)*ThetaScalingFactor+PI;  //--unit: radian
                cell::GeneratePhi(&phi);
                cell::CalculateVelFromThetaPhi(speed, &theta, &phi, Vel);
                *energy = 0.5* *mass*(*speed)*(*speed); //--unit: Joule
        }else if (*ParticleType == iClIonType){
                double ThetaScalingFactor = 1;
                cell::GenerateIonThetaEnergy(cumulativeflux_ClIon, IonAngle, IonEnergy, &theta, energy);
                theta = (theta-PI)*IonThetaScalingFactor+PI; //--unit: radian
                *speed = sqrt(2* *energy/(*mass));
                cell::GeneratePhi(&phi);
                cell::CalculateVelFromThetaPhi(speed, &theta, &phi, Vel);
        }else if (*ParticleType == iCl2IonType){
                double ThetaScalingFactor = 1;
                cell::GenerateIonThetaEnergy(cumulativeflux_Cl2Ion, IonAngle, IonEnergy, &theta, energy);
                theta = (theta - PI)*IonThetaScalingFactor+PI; //--unit: radian
                *speed = sqrt(2* *energy/(*mass));
                cell::GeneratePhi(&phi);
                cell::CalculateVelFromThetaPhi(speed, &theta, &phi, Vel);
        }else if (*ParticleType == iArIonType){
                double ThetaScalingFactor = 1;
                cell::GenerateIonThetaEnergy(cumulativeflux_ArIon, IonAngle, IonEnergy, &theta, energy);
                theta = (theta-PI)*IonThetaScalingFactor+PI; //--unit: radian
                *speed = sqrt(2**energy/(*mass));
                cell::GeneratePhi(&phi);
                cell::CalculateVelFromThetaPhi(speed, &theta, &phi, Vel);
        }
}

//=============================================================================
void cell::Propagation(double* Temperature, int* PropagationTimestepNumber, double* ParticleSizeFactor,
                                                int SurfaceSearchingIndex [][3], int* SurfaceSearchingNumber, int* ParticleType, double* mass, double* dPos, int* iPos,
                                               double* Vel, double* speed, double* energy)
{
        int iPreviousTag;                                                                         //--Previous iTag of Particle Center
        double dPreviousPos [3];                                                       //--Previous dPos (x, y, z) of Particle Center
        int iPreviousPos [3];                                                                  //--Previous iPos (x, y, z) of Particle Center

        int iCurrentTag;                                                                           //--Current iTag of Particle Center
        int iCurrentEdgeTag [6] ;                                                        //--Current iTag of Particle's neighboring six Edge points
        int iCurrentEdgePos [6][3];                                                   //--Current iPos (x, y, z) of Particle's neighboring six Edge points
        int CountPointInSolid = 0;                                                     //--to count how many point of a seven-point molecule is on solid cell


        //cout << "iPos =" << iPos[X_dir] << " " << iPos[Y_dir] << " " << iPos[Z_dir] << endl;
        //cout << "During Propagation : "<< endl;

        int ParticleInDomain = 1;
        int iCollisionTag = 0;
        int iPreviousCollisionTag = 0;
        double PropagationTimeInterval = cell::dDimLength_per_cell[X_dir]/(*speed);  //--Time Interval for propagation ;

        for(int indexStep = 0; indexStep < *PropagationTimestepNumber; indexStep++){
                if( *ParticleType == 0){
                        break;
                }

                iPreviousPos[X_dir] = iPos[X_dir];
                iPreviousPos[Y_dir] = iPos[Y_dir];
                iPreviousPos[Z_dir] = iPos[Z_dir];
                dPreviousPos[X_dir] =dPos[X_dir] ;
                dPreviousPos[Y_dir] =dPos[Y_dir] ;
                dPreviousPos[Z_dir] =dPos[Z_dir] ;
                iPreviousTag = iPreviousPos[X_dir] + ( iPreviousPos[Y_dir] + iPreviousPos[Z_dir]*cell::iDimSize[Y_dir] )*cell::iDimSize[X_dir];

                //--Propagation
                dPos[X_dir]   += Vel[X_dir]*PropagationTimeInterval ;
                dPos[Y_dir]   += Vel[Y_dir]*PropagationTimeInterval;
                dPos[Z_dir]   += Vel[Z_dir]*PropagationTimeInterval;
                iPos[X_dir]     = floor(  dPos[X_dir]/cell::dDimLength_per_cell[X_dir]);
                iPos[Y_dir]     = floor(  dPos[Y_dir]/cell::dDimLength_per_cell[Y_dir]);
                iPos[Z_dir]     = floor(  dPos[Z_dir]/cell::dDimLength_per_cell[Z_dir]);
                cell::ParticleEdge(ParticleSizeFactor, dPos, iCurrentEdgePos, &ParticleInDomain );
                cell::BoundaryMapping(iPos, dPos, &ParticleInDomain);

                if(ParticleInDomain == 0){
                        break;
                }


                //--Count how many point for a seven-point molecule is on solid
                CountPointInSolid = 0;
                for(int i=0; i<6; i++){
                        iCurrentEdgeTag[i] =  iCurrentEdgePos[i][X_dir]
                                                                         +( iCurrentEdgePos[i][Y_dir] + iCurrentEdgePos[i][Z_dir]*cell::iDimSize[Y_dir] )*cell::iDimSize[X_dir];
                        if (cell::iStatus[iCurrentEdgeTag[i]] == iSubstrateStat || cell::iStatus[iCurrentEdgeTag[i]] == iMaskStat)    CountPointInSolid++;
                }
                iCurrentTag =  iPos[X_dir] + ( iPos[Y_dir] + iPos[Z_dir]*iDimSize[Y_dir] )*iDimSize[X_dir];
                if (cell::iStatus[iCurrentTag]  == iSubstrateStat || cell::iStatus[iCurrentTag] == iMaskStat)    CountPointInSolid++;
                //cout << "CountPointInSolid = " << CountPointInSolid << endl;

                //--Check if a particle collide on solid cell
                if (CountPointInSolid < 1){
                        if(Vel[Z_dir] == 0){
                                break;
                        }
                        continue;
                }else if (CountPointInSolid  > 1){
                        dPos[X_dir] = dPreviousPos[X_dir];
                        dPos[Y_dir] = dPreviousPos[Y_dir];
                        dPos[Z_dir] = dPreviousPos[Z_dir];
                        PropagationTimeInterval = 0.5 * PropagationTimeInterval;
                        if ( PropagationTimeInterval*(*speed)/cell::dDimLength_per_cell[X_dir] < 1E-10){
                                break;
                        }
                        continue;
                }else if (CountPointInSolid == 1){
                        for (int i = 0; i < 6 ; i++){
                                if ( iStatus[iCurrentEdgeTag[i]] == iSubstrateStat || iStatus[iCurrentEdgeTag[i]] == iMaskStat){
                                        iPreviousCollisionTag = iCurrentTag;
                                        iCollisionTag = iCurrentEdgeTag[i];
                                         break;
                                }
                        }
                        //--To record the normalized normal vector and reflected vector
                        double normSurfaceNormal [3] = {0.0};               //--normalized surface normal vector (x, y, z)
                        double normIncidentVelocity [3] = {0.0};
                        double normReflectedVelocity [3] = {0.0};          //--normalized reflected velocity vector (x, y, z)
                        double GrazingAngle = 0.0;                                         //--angle between surface and velocity
                        double IncidentAngle = 0.0;                                         //--angle between normal and velocity
                        double NdotV;
                        double arccos_NdotV;
                        //--Output: normSurfaceNormal, normReflectedVelocity, GrazingAngle, IncidentAngle
                        cell::SurfaceNormal(SurfaceSearchingIndex, SurfaceSearchingNumber, iPos, Vel, speed, normSurfaceNormal);
                        cell::normReflectedVelocity(normSurfaceNormal, Vel, speed, &NdotV, normReflectedVelocity);
                        cell::SurfaceAngle(&NdotV, ParticleType, &IncidentAngle, &GrazingAngle);
                        if(*ParticleType == 0){
                                break;
                        }

                        if ( cell::iStatus[iCollisionTag] == iMaskStat){
                                //--Output: ParticleType, mass
                                cell::MaskReaction( &iCollisionTag, energy, &IncidentAngle, ParticleType, mass);
                                if(*ParticleType == 0){
                                        break;
                                }
                                if ( *ParticleType == iClIonType || *ParticleType == iCl2IonType || *ParticleType == iArIonType)    break;
                                cell::MaskReflection( normReflectedVelocity, ParticleType, Vel, speed );
                                if ( *speed == 0  ){
                                        break;
                                }
                        }else if ( cell::iStatus[iCollisionTag] == iSubstrateStat){
                                //--Output: ParticleType, mass
                                cell::SubstrateReaction( &iCollisionTag, &iPreviousCollisionTag, energy, &IncidentAngle, ParticleType, mass);
                                if( *ParticleType == 0){
                                        break;
                                }

                                //--Output: Vel, speed, energy
                                cell::SubstrateReflection( Temperature, normSurfaceNormal, normReflectedVelocity, &GrazingAngle,
                                                                                    ParticleType, mass, Vel, speed, energy );
                                if ( *speed == 0  ){
                                        break;
                                }
                        }
                        continue;
                        //cout << "iParticle = " << iParticle << "  " << "ParticleType = " << ParticleType << endl;
                }
        }//--End of while loop if CountPointInSolid == 1, that is, particle collide with solid (substrate or mask)

}

//=====================================================================================
//--Output: normSurfaceNormal, normReflecteVelocity, GgrazinAangle, IncidenyAngle
void cell::SurfaceNormal(int SearchingIndex [][3], int* SearchingNumber,  int* iPos, double* Vel, double* speed, double* normSurfaceNormal){

        double Surfacesites [*SearchingNumber][3];
        int indexSurfacesites = 0;
        cell::SurfaceSites( SearchingIndex, SearchingNumber, iPos, Surfacesites, &indexSurfacesites  );

        double xbar, ybar, zbar, xsum, ysum, zsum, numSites, A11, A12, A13, A21, A22, A23, A31, A32, A33;
        xbar = 0; ybar = 0; zbar = 0; xsum = 0; ysum = 0; zsum = 0; numSites = 0;
        A11 = 0; A12 = 0; A13 = 0; A21 = 0; A22 = 0; A23 = 0; A31 = 0; A32 = 0; A33 = 0;

        for(int i =0; i< indexSurfacesites; i++){
                xsum += Surfacesites[i][X_dir];
                ysum += Surfacesites[i][Y_dir];
                zsum += Surfacesites[i][Z_dir];
                numSites++;
        }

        xbar = xsum/numSites;
        ybar = ysum/numSites;
        zbar = zsum/numSites;

        //--calculated surface normal
        for(int i =0; i< indexSurfacesites; i++){
                A11 = A11 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][X_dir] - xbar);
                A22 = A22 + (Surfacesites[i][Y_dir] - ybar)*(Surfacesites[i][Y_dir] - ybar);
                A33 = A33 + (Surfacesites[i][Z_dir] - zbar)*(Surfacesites[i][Z_dir] - zbar);
                A12 = A12 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][Y_dir] - ybar);
                A13 = A13 + (Surfacesites[i][X_dir] - xbar)*(Surfacesites[i][Z_dir] - zbar);
                A23 = A23 + (Surfacesites[i][Y_dir] - ybar)*(Surfacesites[i][Z_dir] - zbar);
        }

        A21 = A12;
        A31 = A13;
        A32 = A23;

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
        normSurfaceNormal[X_dir] = V[minVal_index*3+X_dir];
        normSurfaceNormal[Y_dir] = V[minVal_index*3+Y_dir];
        normSurfaceNormal[Z_dir] = V[minVal_index*3+Z_dir];

        //--determine the direction of the surface normal
        int iPositivePoll = 0;
        int iNegativePoll = 0;
        int iPosPoll [3];
        int iTagPoll;

        //cout << "position = " << iPos[X_dir] << " " << iPos[Y_dir] << " " << iPos[Z_dir] << endl;

        for (int i = 1; i <= 3 ; i++){
                iPosPoll[X_dir] = floor(iPos[X_dir] + 2*i*normSurfaceNormal[X_dir]);
                iPosPoll[Y_dir] = floor(iPos[Y_dir] + 2*i*normSurfaceNormal[Y_dir]);
                iPosPoll[Z_dir] = floor(iPos[Z_dir] + 2*i*normSurfaceNormal[Z_dir]);
                //--periodic boundary condition
                while ( iPosPoll[X_dir] >= cell::iDimSize[X_dir] )     iPosPoll[X_dir] -= cell::iDimSize[X_dir];
                while ( iPosPoll[Y_dir] >= cell::iDimSize[Y_dir] )     iPosPoll[Y_dir] -= cell::iDimSize[Y_dir];
                while ( iPosPoll[Z_dir] >= cell::iDimSize[Z_dir] )     iPosPoll[Z_dir] -= cell::iDimSize[Z_dir];
                while ( iPosPoll[X_dir] < 0       )                                            iPosPoll[X_dir] += cell::iDimSize[X_dir];
                while ( iPosPoll[Y_dir] < 0       )                                            iPosPoll[Y_dir] += cell::iDimSize[Y_dir];
                while ( iPosPoll[Z_dir] < 0       )                                            iPosPoll[Z_dir] += cell::iDimSize[Z_dir];
                iTagPoll = iPosPoll[X_dir] + (iPosPoll[Y_dir] + iPosPoll[Z_dir]* cell::iDimSize[Y_dir])*cell::iDimSize[X_dir];
                if ( cell::iStatus[iTagPoll] == iVacuumStat ){
                        iPositivePoll++;
                }
        }

        for (int i = 1; i <= 3 ; i++){
                iPosPoll[X_dir] = floor(iPos[X_dir] - 2*i*normSurfaceNormal[X_dir]);
                iPosPoll[Y_dir] = floor(iPos[Y_dir] - 2*i*normSurfaceNormal[Y_dir]);
                iPosPoll[Z_dir] = floor(iPos[Z_dir] - 2*i*normSurfaceNormal[Z_dir]);
                //--periodic boundary condition
                while ( iPosPoll[X_dir] >= cell::iDimSize[X_dir] )     iPosPoll[X_dir] -= cell::iDimSize[X_dir];
                while ( iPosPoll[Y_dir] >= cell::iDimSize[Y_dir] )     iPosPoll[Y_dir] -= cell::iDimSize[Y_dir];
                while ( iPosPoll[Z_dir] >= cell::iDimSize[Z_dir] )     iPosPoll[Z_dir] -= cell::iDimSize[Z_dir];
                while ( iPosPoll[X_dir] < 0       )                                             iPosPoll[X_dir] += cell::iDimSize[X_dir];
                while ( iPosPoll[Y_dir] < 0       )                                             iPosPoll[Y_dir] += cell::iDimSize[Y_dir];
                while ( iPosPoll[Z_dir] < 0       )                                             iPosPoll[Z_dir] += cell::iDimSize[Z_dir];
                iTagPoll = iPosPoll[X_dir] + (iPosPoll[Y_dir] + iPosPoll[Z_dir]*cell::iDimSize[Y_dir])*cell::iDimSize[X_dir];
                if ( cell::iStatus[iTagPoll] == iVacuumStat ){
                        iNegativePoll++;
                }
        }


        if ( iPositivePoll > iNegativePoll){
                 normSurfaceNormal[X_dir] = (+1)*normSurfaceNormal[X_dir];
                 normSurfaceNormal[Y_dir] = (+1)*normSurfaceNormal[Y_dir];
                 normSurfaceNormal[Z_dir] = (+1)*normSurfaceNormal[Z_dir];
        }else if ( iPositivePoll < iNegativePoll){
                 normSurfaceNormal[X_dir] = (-1)*normSurfaceNormal[X_dir];
                 normSurfaceNormal[Y_dir] = (-1)*normSurfaceNormal[Y_dir];
                 normSurfaceNormal[Z_dir] = (-1)*normSurfaceNormal[Z_dir];
         }else if ( iPositivePoll == iNegativePoll){
                 double randPoll = ran3(&iRandTag);
                 //cout << "rdd in surface normal = " << rdd << endl;
                 if( randPoll < 0.5 ){
                         normSurfaceNormal[X_dir] = (+1)*normSurfaceNormal[X_dir];
                         normSurfaceNormal[Y_dir] = (+1)*normSurfaceNormal[Y_dir];
                         normSurfaceNormal[Z_dir] = (+1)*normSurfaceNormal[Z_dir];
                 }else{
                         normSurfaceNormal[X_dir] = (-1)*normSurfaceNormal[X_dir];
                         normSurfaceNormal[Y_dir] = (-1)*normSurfaceNormal[Y_dir];
                         normSurfaceNormal[Z_dir] = (-1)*normSurfaceNormal[Z_dir];
                  }
         }

}


//--calculate the normalized reflected velocity
void cell::normReflectedVelocity(double* normN, double* Vel, double* speed, double* NdotV, double* normReflectedVelocity){
        double normV [3];
        normV[X_dir] = Vel[X_dir]/(*speed);  //--x component of normalized incident velocity vector
        normV[Y_dir] = Vel[Y_dir]/(*speed);  //--y component of normalized incident velocity vector
        normV[Z_dir] = Vel[Z_dir]/(*speed);  //--z component of normalized incident velocity vector
        *NdotV= normN[X_dir]*normV[X_dir]+normN[Y_dir]*normV[Y_dir]+normN[Z_dir]*normV[Z_dir];
        normReflectedVelocity[X_dir] = normV[X_dir] - 2*(*NdotV)*normN[X_dir];
        normReflectedVelocity[Y_dir] = normV[Y_dir] - 2*(*NdotV)*normN[Y_dir];
        normReflectedVelocity[Z_dir] = normV[Z_dir] - 2*(*NdotV)*normN[Z_dir];
}

void cell::SurfaceAngle(double* NdotV, int* ParticleType, double* IncidentAngle, double* GrazingAngle){
        double arccos_NdotV = acos(*NdotV)*180/PI;
        if ( arccos_NdotV  <  90 ){
                *ParticleType = 0;
        }
        *IncidentAngle = 180 - arccos_NdotV;
        *GrazingAngle =  90 - *IncidentAngle;
}

//==========================================================================
void cell::MaskReaction( int* iCollisionTag, double* energy, double* IncidentAngle, int* ParticleType, double* mass)
{

        int ReactionIndex;
        int EmittedParticle = 0;                                   //--index for emitted particle such as SiClx(g)
        int ReflectedParticle= 0;                              //--index for original reflected particle such as Ar*, Cl*

        if (  *ParticleType == iClRadicalType ){        ReflectedParticle = iClRadicalType;
        }else if (  *ParticleType == iSiClgType ){     ReflectedParticle = iSiClgType;
        }else if (  *ParticleType == iSiCl2gType ){  ReflectedParticle = iSiCl2gType;
        }else if (  *ParticleType == iSiCl3gType ){  ReflectedParticle = iSiCl3gType;
        }else if (  *ParticleType == iClIonType ){     ReflectedParticle = iClIonType;
        }else if (  *ParticleType == iCl2IonType ){  ReflectedParticle = iCl2IonType;
        }else if (   *ParticleType == iArIonType ){   ReflectedParticle = iArIonType;
        }

        //--Assinge Product Particle Type and mass
        *ParticleType = ReflectedParticle;
        if( *ParticleType == iClRadicalType){           *mass = MassChlorine;
        }else if ( *ParticleType == iSiClgType){       *mass = MassSilicon+MassChlorine;
        }else if ( *ParticleType == iSiCl2gType){    *mass = MassSilicon+2*MassChlorine;
        }else if ( *ParticleType == iSiCl3gType){    *mass = MassSilicon+3*MassChlorine;
        }else if (  *ParticleType == iClIonType){      *mass = MassChlorine;
        }else if ( *ParticleType == iCl2IonType){    *mass = MassChlorine*2;
        }else if (  *ParticleType == iArIonType){     *mass = MassArgon;
        }

}



void cell::SubstrateReaction( int* iCollisionTag, int* iPreviousCollisionTag, double* energy, double* IncidentAngle, int* ParticleType, double* mass){
        int ReactionIndex;
        int EmittedParticle = 0;                                   //--index for emitted particle such as SiClx(g)
        int ReflectedParticle= 0;                              //--index for original reflected particle such as Ar*, Cl*
        double p0_ClRadicalReaction [6] = {0.99, 0.40, 0.30, 0.02, 0.0001, 0.08};
        double p0_redeposition [3] = {0.02, 0.02, 0.02};

        double p0_ClIonReaction [5] = {0.05, 0.10, 0.20, 0.50, 0.50};
        double Eth_ClIonReaction [5] = {25, 35, 10, 10, 10};
        double E0_ClIonReaction = 100;
        int type_ClIonReaction [5] = {P_sputtering, P_sputtering, C_sputtering, C_sputtering, C_sputtering};

        double p0_Cl2IonReaction [6] = {0.02, 0.20, 0.25, 0.25, 0.25, 0.25};
        double Eth_Cl2IonReaction [6] = {25, 10, 10, 10, 10, 10};
        double E0_Cl2IonReaction = 100;
        int type_Cl2IonReaction [6] = {P_sputtering, C_sputtering, C_sputtering, C_sputtering, C_sputtering, C_sputtering};

        double p0_ArIonReaction [4] = {0.05, 0.20, 0.50, 0.50};
        double Eth_ArIonReaction [4] = {25, 10, 10, 10};
        double E0_ArIonReaction = 100;
        int type_ArIonReaction [4] = {P_sputtering, C_sputtering, C_sputtering, C_sputtering};

        double phys_sputter_prob [10] = {0.55, 0.555, 0.60, 0.555, 0.85, 1.3, 1.355, 1.0, 0.75, 0.0};
        double chem_sputter_prob [10] = {1.00, 1.000, 1.00, 1.000, 1.00, 0.9, 0.800, 0.6, 0.30, 0.0};

        if( *ParticleType == iClRadicalType ){
                int number_of_reactions = 6;
                double ReactionProb [number_of_reactions];


                cell::RadicalReactionProb(ParticleType, iCollisionTag, &number_of_reactions, p0_ClRadicalReaction, ReactionProb);
                cell::ClRadicalReaction(iCollisionTag, ReactionProb, &ReflectedParticle, &ReactionIndex);

        }else if (  *ParticleType == iSiClgType){
                cell::Redeposition(ParticleType, iPreviousCollisionTag, p0_redeposition, &ReflectedParticle, &ReactionIndex);

        }else if (  *ParticleType == iSiCl2gType){
                cell::Redeposition(ParticleType, iPreviousCollisionTag, p0_redeposition, &ReflectedParticle, &ReactionIndex);

        }else if ( *ParticleType == iSiCl3gType){
                cell::Redeposition(ParticleType, iPreviousCollisionTag, p0_redeposition, &ReflectedParticle, &ReactionIndex);

        }else if (  *ParticleType == iClIonType){
                int number_of_reactions = 5;
                double ReactionProb [number_of_reactions];
                cell::IonReactionProb(ParticleType, iCollisionTag,  &number_of_reactions, Eth_ClIonReaction, &E0_ClIonReaction, p0_ClIonReaction,
                                                                type_ClIonReaction, phys_sputter_prob, chem_sputter_prob, energy, IncidentAngle, ReactionProb );
                cell::ClIonReaction(iCollisionTag, energy, IncidentAngle, ReactionProb, &number_of_reactions,
                                                        &ReflectedParticle, &EmittedParticle, &ReactionIndex);
        }else if (  *ParticleType == iCl2IonType){
                int number_of_reactions = 6;
                double ReactionProb [number_of_reactions];
                cell::IonReactionProb(ParticleType, iCollisionTag,  &number_of_reactions, Eth_Cl2IonReaction, &E0_Cl2IonReaction, p0_Cl2IonReaction,
                                                                type_Cl2IonReaction, phys_sputter_prob, chem_sputter_prob, energy, IncidentAngle, ReactionProb );
                cell::Cl2IonReaction(iCollisionTag, energy, IncidentAngle, ReactionProb, &number_of_reactions,
                                                            &ReflectedParticle, &EmittedParticle, &ReactionIndex);
        }else if (  *ParticleType == iArIonType){
                int number_of_reactions = 4;
                double ReactionProb [number_of_reactions];
                cell::IonReactionProb(ParticleType, iCollisionTag,  &number_of_reactions, Eth_ArIonReaction, &E0_ArIonReaction, p0_ArIonReaction,
                                                                type_ArIonReaction, phys_sputter_prob, chem_sputter_prob, energy, IncidentAngle, ReactionProb );
                cell::ArIonReaction(iCollisionTag, energy, IncidentAngle, ReactionProb, &number_of_reactions,
                                                        &ReflectedParticle, &EmittedParticle, &ReactionIndex);

        }

        //--Assinge Product Particle Type and mass
        *ParticleType = ReflectedParticle;
        if( *ParticleType == iClRadicalType){           *mass = MassChlorine;
        }else if ( *ParticleType == iSiClgType){       *mass = MassSilicon+MassChlorine;
        }else if ( *ParticleType == iSiCl2gType){    *mass = MassSilicon+2*MassChlorine;
        }else if ( *ParticleType == iSiCl3gType){    *mass = MassSilicon+3*MassChlorine;
        }else if (  *ParticleType == iClIonType){      *mass = MassChlorine;
        }else if ( *ParticleType == iCl2IonType){    *mass = MassChlorine*2;
        }else if (  *ParticleType == iArIonType){     *mass = MassArgon;
        }
}




void cell::MaskReflection( double* normReflectedVelocity, int* ParticleType, double* Vel, double* speed )
{
        if( *ParticleType == iClRadicalType ){//--iClRadical, iSiClg, iSiCl2g, iSiCl3g
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iSiClgType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iSiCl2gType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iSiCl3gType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iClIonType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iCl2IonType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }else if ( *ParticleType == iArIonType ){
                cell::ElasticReflection(normReflectedVelocity, speed, Vel );
        }
}




void cell::SubstrateReflection( double* Temperature, double* normSurfaceNormal, double* normReflectedVelocity,
                                                             double* GrazingAngle, int* ParticleType, double* mass, double* Vel, double* speed, double* energy )
{
        if( *ParticleType == iClRadicalType ){//--iClRadical, iSiClg, iSiCl2g, iSiCl3g
                double ReemissionCosineLawPower = 100;
                cell::DiffusiveReflection(Temperature, normSurfaceNormal, &ReemissionCosineLawPower, mass, Vel, speed, energy);
        }else if ( *ParticleType == iSiClgType ){
                double ReemissionCosineLawPower = 100;
                cell::DiffusiveReflection(Temperature, normSurfaceNormal, &ReemissionCosineLawPower, mass, Vel, speed, energy);
        }else if ( *ParticleType == iSiCl2gType ){
                double ReemissionCosineLawPower = 100;
                cell::DiffusiveReflection(Temperature, normSurfaceNormal, &ReemissionCosineLawPower, mass, Vel, speed, energy);
        }else if ( *ParticleType == iSiCl3gType ){
                double ReemissionCosineLawPower = 100;
                cell::DiffusiveReflection(Temperature, normSurfaceNormal, &ReemissionCosineLawPower, mass, Vel, speed, energy);
        }else if ( *ParticleType == iClIonType ){
                cell::InelasticReflection(normReflectedVelocity, GrazingAngle, mass, Vel, speed, energy);
        }else if ( *ParticleType == iCl2IonType ){
                cell::InelasticReflection(normReflectedVelocity, GrazingAngle, mass, Vel, speed, energy);
        }else if ( *ParticleType == iArIonType ){
                cell::InelasticReflection(normReflectedVelocity, GrazingAngle, mass, Vel, speed, energy);
        }
}
