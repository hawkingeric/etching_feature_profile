#include "etching.h"
#include "cell.h"
#include "rand.h"

using namespace std;
//--Output: ParticleType and mass, Input: CumuProb
void cell::GenerateParticleTypeMass(double* CumuProb, int* ParticleType, double* mass){
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randParticle = unif(generator); //--Random number for Particle Type
        if (randParticle >= 0 && randParticle < CumuProb[0]){
                *ParticleType = iClRadicalType;
                *mass = MassChlorine;
        }else if ( randParticle >= CumuProb[0] && randParticle < CumuProb[1]){
                *ParticleType = iClIonType;
                *mass = MassChlorine;
        }else if ( randParticle >= CumuProb[1] && randParticle < CumuProb[2]){
                *ParticleType = iCl2IonType;
                *mass = MassChlorine*2;
        }else if ( randParticle >= CumuProb[2] && randParticle < CumuProb[3]){
                *ParticleType = iArIonType;
                *mass = MassArgon;
        }else{
                *ParticleType = 0;
                cout << "ParticleType = " << *ParticleType << endl;
        }
}

//--Output: dPos and iPos
void cell::GeneratePosition( double* dPos, int* iPos){
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randPosX = unif(generator);  //--Random number for X coordinate
        double randPosY = unif(generator);  //--Random number for Y coordinate
        dPos[X_dir] = randPosX*cell::dDimLength[X_dir];
        dPos[Y_dir] = randPosY*cell::dDimLength[Y_dir];
        dPos[Z_dir] = cell::dDimLength[Z_dir]-cell::dDimLength_per_cell[Z_dir]*0.5;
        //cout << "Nx Ny Nz = " << cell::iDimSize[X_dir] << " " << cell::iDimSize[Y_dir] << " " << cell::iDimSize[Z_dir] << endl;
        //cout << "dx dy dz = " << cell::dDimLength_per_cell[X_dir] << " " << cell::dDimLength_per_cell[Y_dir] << " " << cell::dDimLength_per_cell[Z_dir] << endl;
        //cout << "Lx Ly Lz = " << cell::dDimLength[X_dir] << " " << cell::dDimLength[Y_dir] << " " << cell::dDimLength[Z_dir] << endl;
        iPos[X_dir]     = floor(  dPos[X_dir]/cell::dDimLength_per_cell[X_dir]);
        iPos[Y_dir]     = floor(  dPos[Y_dir]/cell::dDimLength_per_cell[Y_dir]);
        //cout << " dPos[Z_dir]/cell::dDimLength_per_cell[Z_dir] = " << dPos[Z_dir]/cell::dDimLength_per_cell[Z_dir] << endl;
        iPos[Z_dir]     = floor(  dPos[Z_dir]/cell::dDimLength_per_cell[Z_dir]);
        //cout << iPos[Z_dir] << endl;

}

//--Output: Radical Velocity and speed
void cell::GenerateRadicalSpeed(double* Temperature, double* v_cut, double* mass, double* speed){
        random_device rd;
        std::default_random_engine generator( rd() ) ;
        uniform_real_distribution<double> unif(0.0, 1.0);  // --Random number uniform distribution
        double v_mp = sqrt(2*BoltzmannConstant_SI* (*Temperature)/(*mass) );  //--the most probable velocity = sqrt(2kT/m)
        double p_max = Maxwell_Boltzmann_pdf(*mass, *Temperature, v_mp);  //--the maximal probability = p(v_mp)
        double v_rand, p_rand, prob_density;
        do
        {
                v_rand = unif(generator)*(*v_cut);          //--Generate first random numbr/
                prob_density = Maxwell_Boltzmann_pdf(*mass, *Temperature, v_rand);   //--Mawwell Boltzmann distribution probability density funcion
                p_rand = unif(generator)*p_max;         //--Generate second random number
        }while (prob_density < p_rand);
        *speed = v_rand;
}


//--Probability density function of Maxwell Boltzmann distribution
double cell::Maxwell_Boltzmann_pdf(double mass, double T, double v){
        double v_p = sqrt(2*BoltzmannConstant_SI*T/mass); // the most probable speed
        return 4*PI*v*v/pow(PI, 1.5)/pow(v_p, 3.0)*exp( -pow(v/v_p, 2.0)  );
}


//--Cumulative density function of Maxwell Boltzmann distribution
double cell::Maxwell_Boltzmann_cdf(double mass, double T, double v){
        double v_p = sqrt(2*BoltzmannConstant_SI*T/mass); // the most probable speed
        return erf(v/v_p) - 2/sqrt(PI)*(v/v_p)*exp( -pow(v/v_p, 2.0)  );
}


//--Calculate velocity of Cumulative density function of Maxwell Boltzmann distribution
double cell::CalculateSpeedCutoff(double mass, double T){
        double criterion = 0.0001;
        double dv = 0.1;
        double v_rms = sqrt(3*BoltzmannConstant_SI*T/mass);  // the root mean square velocity = sqrt(3kT/m)
        double v = v_rms;
        while(  (1 - Maxwell_Boltzmann_cdf(mass, T, v)  ) > criterion){
                v = v + dv;
        };
        return v;
}




//--Output: Theta, Input: Velocity
void cell::GenerateTheta(double* theta){
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        *theta = unif(generator)*PI/2+PI/2;
}

//--Output: Phi
void cell::GeneratePhi(double* phi){
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randPhi = unif(generator);  //--Random number for Phi
        *phi = randPhi*2*PI;
}

//--Output: Theta and Energy
void cell::GenerateIonThetaEnergy(vector<double>& cumulativeflux, vector<double>& IonAngle,  vector<double>& IonEnergy,
                                                                                double* theta, double* energy)
{
        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double rand1 = unif(generator);  //--Random number for flux
        double rand2 = unif(generator);  //--Random number for intropolation
        for(int i = 0;    i < cumulativeflux.size()-1;     i++){
                if  (  rand1 >= cumulativeflux[i] && rand1 < cumulativeflux[i+1]  ){
                        *theta =(    IonAngle[i] + rand2*(  IonAngle[i+1] - IonAngle[i]  )    )*PI/180+PI;  //--unit: radian
                        *energy = (    IonEnergy[i] + rand2*(  IonEnergy[i+1] - IonEnergy[i]  )    )/Joule_to_eV;  //--unit: Joule
                        break;
                }
        }
}

//--Output: Velocity, Input: Speed, Theta, and Phi
void cell::CalculateVelFromThetaPhi(double* speed, double* theta, double* phi, double* Vel){
        Vel[X_dir] = *speed*sin(*theta)*cos(*phi);
        Vel[Y_dir] = *speed*sin(*theta)*sin(*phi);
        Vel[Z_dir] = *speed*cos(*theta);
}








//====================================================================================//

//--Output: iPos and dPos. For periodic boundary condition, velocities need not change
void cell::BoundaryMapping(int* iPos, double* dPos, int* ParticleInDomain){

        //--Periodic boundary condition
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
        if(  iPos[Z_dir] < 0 || iPos[Z_dir] >= cell::iDimSize[Z_dir] ){
                *ParticleInDomain = 0;
        }
}


//--Output: iParticleEdgePos
void cell::ParticleEdge(double* ParticleSizeFactor, double* dPos, int iParticleEdgePos [][3], int* ParticleInDomain )
{
        //--Pre-calculation of particle size and initialization of six point particle
        double indexParticleEdge [6][3] = {  {-1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, +1, 0}, {0, 0, -1}, {0, 0, +1} };
        double dParticleEdgePos [3];                              //--dPos (x, y, z) of particle's neighboring edge points


        //--Calculate particle neighboring six points (particle size edge)
        for ( int i = 0; i < 6 ; i++){
                dParticleEdgePos[X_dir]   = dPos[X_dir] + indexParticleEdge[i][X_dir]*(*ParticleSizeFactor)*cell::dDimLength_per_cell[X_dir];;
                dParticleEdgePos[Y_dir]   = dPos[Y_dir] + indexParticleEdge[i][Y_dir]*(*ParticleSizeFactor)*cell::dDimLength_per_cell[Y_dir];
                dParticleEdgePos[Z_dir]   = dPos[Z_dir] + indexParticleEdge[i][Z_dir]*(*ParticleSizeFactor)*cell::dDimLength_per_cell[Z_dir];
                iParticleEdgePos[i][X_dir]     = floor(dParticleEdgePos[X_dir]/cell::dDimLength_per_cell[X_dir]  );
                iParticleEdgePos[i][Y_dir]     = floor(dParticleEdgePos[Y_dir]/cell::dDimLength_per_cell[Y_dir]  );
                iParticleEdgePos[i][Z_dir]     = floor(dParticleEdgePos[Z_dir]/cell::dDimLength_per_cell[Z_dir] );

                if( iParticleEdgePos[i][X_dir] < 0){
                        iParticleEdgePos[i][X_dir] += cell::iDimSize[X_dir];
                }else if(  iParticleEdgePos[i][X_dir] >= cell::iDimSize[X_dir]){
                        iParticleEdgePos[i][X_dir] -= cell::iDimSize[X_dir];
                }
                if( iParticleEdgePos[i][Y_dir] < 0){
                        iParticleEdgePos[i][Y_dir] += cell::iDimSize[Y_dir];
                }else if( iParticleEdgePos[i][Y_dir] >= cell::iDimSize[Y_dir]){
                        iParticleEdgePos[i][Y_dir] -= cell::iDimSize[Y_dir];
                }
                if( iParticleEdgePos[i][Z_dir] < 0 || iParticleEdgePos[i][Z_dir] >= cell::iDimSize[Z_dir]){
                        *ParticleInDomain = 0;
                }
        }
}




//=================================================================================//

void cell::SurfaceSites(int SearchingIndex [][3] , int* SearchingNumber, int* iPos, double Surfacesites [][3], int* indexSurfacesites ){
        int iPosNB [3];
        int iPosNBTag [3];
        int  iNBTag;
        *indexSurfacesites = 0;

        //--determine surface sites
        for (int i = 0; i < *SearchingNumber; i++){
                iPosNB[X_dir] = iPos[X_dir] + SearchingIndex[i][X_dir];
                iPosNB[Y_dir] = iPos[Y_dir] + SearchingIndex[i][Y_dir];
                iPosNB[Z_dir] = iPos[Z_dir] + SearchingIndex[i][Z_dir];
                iPosNBTag[X_dir] = iPosNB[X_dir];
                iPosNBTag[Y_dir] = iPosNB[Y_dir];
                iPosNBTag[Z_dir] = iPosNB[Z_dir];

                //--periodic boundary condition
                while (iPosNBTag[X_dir] >= cell::iDimSize[X_dir] )  iPosNBTag[X_dir] -= cell::iDimSize[X_dir];
                while (iPosNBTag[Y_dir] >= cell::iDimSize[Y_dir] )  iPosNBTag[Y_dir] -= cell::iDimSize[Y_dir];
                while (iPosNBTag[Z_dir] >= cell::iDimSize[Z_dir] )  iPosNBTag[Z_dir] -= cell::iDimSize[Z_dir];
                while (iPosNBTag[X_dir] < 0                )                                iPosNBTag[X_dir] += cell::iDimSize[X_dir];
                while (iPosNBTag[Y_dir] < 0                )                                iPosNBTag[Y_dir] += cell::iDimSize[Y_dir];
                while (iPosNBTag[Z_dir] < 0                )                                iPosNBTag[Z_dir] += cell::iDimSize[Z_dir];

                //--calculated itag of Neighboring cell
                iNBTag = iPosNBTag[X_dir] + ( iPosNBTag[Y_dir] + iPosNBTag[Z_dir] * cell::iDimSize[Y_dir])*cell::iDimSize[X_dir];

                //cout << "i = " << i << endl;
                //--check if the cell property of center cell and neighboring cell same
                if ( cell::iStatus[iNBTag] == iSubstrateStat || cell::iStatus[iNBTag] == iMaskStat){
                        //--check six faces of neighboring cells are exposed to vacuum
                        if (  cell::iStatus[cell::iID_NBR[iNBTag][4]] == iVacuumStat || cell::iStatus[cell::iID_NBR[iNBTag][12]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[iNBTag][10]] == iVacuumStat || cell::iStatus[cell::iID_NBR[iNBTag][14]] == iVacuumStat ||
                                cell::iStatus[cell::iID_NBR[iNBTag][16]] == iVacuumStat || cell::iStatus[cell::iID_NBR[iNBTag][22]] == iVacuumStat){
                                Surfacesites[*indexSurfacesites][X_dir] = double(iPosNB[X_dir]);
                                Surfacesites[*indexSurfacesites][Y_dir] = double(iPosNB[Y_dir]);
                                Surfacesites[*indexSurfacesites][Z_dir] = double(iPosNB[Z_dir]);
                                *indexSurfacesites = *indexSurfacesites + 1;
                        }
                }
        }
}



//====================================================================================//

void cell::RadicalReactionProb(int* ParticleType, int* iTag,  int* number_of_reactions, double* p0_ClRadicalReaction, double* ReactionProb )
{
        double sumReactionProb = 0;
        double denominator = cell::dNumSiClxs[*iTag][0] + cell::dNumSiClxs[*iTag][1] + cell::dNumSiClxs[*iTag][2] * 2 + cell::dNumSiClxs[*iTag][3] * 2;
        ReactionProb[0] = p0_ClRadicalReaction[0]*cell::dNumSiClxs[*iTag][0]/denominator;
        ReactionProb[1] = p0_ClRadicalReaction[1]*cell::dNumSiClxs[*iTag][1]/denominator;
        ReactionProb[2] = p0_ClRadicalReaction[2]*cell::dNumSiClxs[*iTag][2]/denominator;
        ReactionProb[3] = p0_ClRadicalReaction[3]*cell::dNumSiClxs[*iTag][2]/denominator;
        ReactionProb[4] = p0_ClRadicalReaction[4]*cell::dNumSiClxs[*iTag][3]/denominator;
        ReactionProb[5] = p0_ClRadicalReaction[5]*cell::dNumSiClxs[*iTag][3]/denominator;
        for( int i = 0; i < *number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }
        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < *number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < *number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }
}

void cell::IonReactionProb(int* ParticleType, int* iTag,  int* number_of_reactions, double* Eth_IonReaction, double* E0,
                                                        double* p0_IonReaction, int* type_IonReaction, double* phys_sputter_prob, double* chem_sputter_prob,
                                                        double* energy, double* IncidentAngle, double* ReactionProb )
{
/*
if(*ParticleType == iClIonType){
        for(int i = 0; i < *number_of_reactions ; i++){
                cout << p0_IonReaction[i] << endl;
                cout << Eth_IonReaction[i] << endl;
                cout << type_IonReaction[i] << endl;
        }
}
*/
        //--Calculation of prob of energy and prob of angle
        double prob_of_energy [*number_of_reactions];
        double prob_of_angle [*number_of_reactions];
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;
        double energy_in_unit_eV = *energy*Joule_to_eV;


        for (int i = 0; i < *number_of_reactions; i++){
                if(energy_in_unit_eV >= Eth_IonReaction[i]){
                        EnhanceFactor = sqrt( (energy_in_unit_eV - Eth_IonReaction[i])/(*E0 - Eth_IonReaction[i]) );
                        prob_of_energy[i] = p0_IonReaction[i]*EnhanceFactor;
                }else if (energy_in_unit_eV < Eth_IonReaction[i]){
                        prob_of_energy[i] = 0.0;
                }

                if( type_IonReaction[i] == P_sputtering ){
                        angle_index = int( floor(*IncidentAngle/10)  );
                        angle_interval = 10.0;
                        angle = double(angle_index)*angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( *IncidentAngle == angle ){
                                prob_of_angle[i] = phys_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - *IncidentAngle)/(*IncidentAngle - angle);
                                prob_of_angle[i] = (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }else if (type_IonReaction[i] == C_sputtering){
                        angle_index = int(floor(*IncidentAngle/10));
                        angle_interval = 10.0;
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  *IncidentAngle == angle ){
                                prob_of_angle[i] = chem_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - *IncidentAngle)/(*IncidentAngle - angle);
                                prob_of_angle[i] = (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }
        }



        double denominator;
        if(*ParticleType == iClIonType){
                denominator = cell::dNumSiClxs[*iTag][0] + cell::dNumSiClxs[*iTag][1]*2 + cell::dNumSiClxs[*iTag][2] + cell::dNumSiClxs[*iTag][3] ;
        }else if (*ParticleType == iCl2IonType){
                denominator = cell::dNumSiClxs[*iTag][0] + cell::dNumSiClxs[*iTag][1] + cell::dNumSiClxs[*iTag][2]*2 + cell::dNumSiClxs[*iTag][3]*2 ;
        }else if (*ParticleType == iArIonType){
                denominator = cell::dNumSiClxs[*iTag][0] + cell::dNumSiClxs[*iTag][1] + cell::dNumSiClxs[*iTag][2] + cell::dNumSiClxs[*iTag][3] ;
        }

        for(int i = 0; i < *number_of_reactions; i++){
                ReactionProb[i] = prob_of_energy[i]*prob_of_angle[i]*cell::dNumSiClxs[*iTag][i]/denominator;
        }

        double sumReactionProb = 0;
        for( int i = 0; i < *number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }

        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < *number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < *number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }
}

void cell::ClRadicalReaction(int* iTag, double* ReactionProb, int* ReflectedParticle, int* ReactionIndex){

        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randReaction = unif(generator);
        if  ( randReaction < ReactionProb[0]  ){
                //--reaction 1 : Si(s) + Cl --> SiCl(s)        p0_ClRadicalReaction[0] = 0.99
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][0]--;
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]++;
                *ReflectedParticle = 0;
                *ReactionIndex = 1;
        }else if  (  randReaction >= ReactionProb[0] && randReaction < ReactionProb[1]  ){
                //--reaction 2 : SiCl(s) + Cl --> SiCl2(s)        p0_ClRadicalReaction[1] = 0.40
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]--;
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]++;
                *ReflectedParticle = 0;
                *ReactionIndex = 2;
        }else if  (  randReaction >= ReactionProb[1] && randReaction < ReactionProb[2]  ){
                //--reaction 3 : SiCl2(s) + Cl --> SiCl3(s)              p0_ClRadicalReaction[2] = 0.30
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]++;
                *ReflectedParticle = 0;
                *ReactionIndex = 3;
        }else if ( randReaction >= ReactionProb[2] && randReaction < ReactionProb[3]   ){
                //--reaction 4 : SiCl2(s) + Cl --> SiCl(s) + Cl2         p0_ClRadicalReaction[3] = 0.02
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]++;
                *ReflectedParticle = 0;
                *ReactionIndex = 4;
        }else if  ( randReaction >= ReactionProb[3] && randReaction < ReactionProb[4]  ){
                //--reaction 5 : SiCl3(s) + Cl --> SiCl4(g)                p0_ClRadicalReaction[4] = 0.0001
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = 0;
                *ReactionIndex = 5;
        }else if ( randReaction >= ReactionProb[4] && randReaction < ReactionProb[5]  )  {
                //--reaction 6 : SiCl3(s) + Cl --> SiCl2(s) + Cl2    p0_ClRadicalReaction[5] = 0.08
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]++;
                *ReflectedParticle = 0;
                *ReactionIndex = 6 ;
        }else{
                *ReflectedParticle = iClRadicalType;
                *ReactionIndex = 0;
        }
        if ( cell::dNumMaterial[*iTag] == 0 )   cell::setStatus(*iTag, iVacuumStat, 1);
}


void cell::Redeposition( int* AdsorbParticle,  int* iTag, double* p0_redeposition, int* ReflectedParticle, int* ReactionIndex){

        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randReaction = unif(generator);
        //--reaction 8 : M(s) + SiClx(g)  --> M(s) + SiClx(s)   p = 0.02   x=1,2,3
        if ( *AdsorbParticle == iSiClgType ){
                if ( randReaction < p0_redeposition[0] ){
                        //--reaction 8-1 : M(s)  +   SiCl(g)  -->  M(s)  +  SiCl(s)        p0_redeposition[0] = 0.02
                        if ( cell::iStatus[*iTag] == iVacuumStat)   cell::iStatus[*iTag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[*iTag][1]++;
                        #pragma omp atomic
                        cell::dNumMaterial[*iTag]++;
                        *ReflectedParticle = 0;
                        *ReactionIndex = 8;
                }else{
                        *ReflectedParticle = iSiClgType;
                }
        }else if ( *AdsorbParticle == iSiCl2gType ){
                if(  randReaction < p0_redeposition[1] ){
                        //--reaction 8-2 : M(s) +   SiCl2(g)  -->  M(s)  +  SiCl2(s)        p0_redeposition[1] = 0.02
                        if ( cell::iStatus[*iTag] == iVacuumStat  )   cell::iStatus[*iTag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[*iTag][2]++;
                        #pragma omp atomic
                        cell::dNumMaterial[*iTag]++;
                        *ReflectedParticle = 0;
                        *ReactionIndex = 8;
                }else{
                        *ReflectedParticle = iSiCl2gType;
                }
        }else if (  *AdsorbParticle == iSiCl3gType ){
                if (  randReaction < p0_redeposition[2] ){
                        //--reaction 8-3 : M(s)  +  SiCl3(g)  -->  M(s)  +  SiCl3(s)        p0_redeposition[2] = 0.02
                        if ( cell::iStatus[*iTag] == iVacuumStat  )   cell::iStatus[*iTag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[*iTag][3]++;
                        #pragma omp atomic
                        cell::dNumMaterial[*iTag]++;
                        *ReflectedParticle = 0;
                        *ReactionIndex = 8;
                }else{
                        *ReflectedParticle = iSiCl3gType;
                }
        }
}



void cell::ClIonReaction(int* iTag, double* energy, double* IncidentAngle, double* ReactionProb, int* number_of_reactions,
                                                    int* ReflectedParticle, int* EmitParticle, int* ReactionIndex){

        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randReaction = unif(generator);
        if ( randReaction < ReactionProb[0]   ){
                //--reaction 9 : Si(s) + Cl+ --> Si(g) + Cl*        p0_ClIonReaction[0] = 0.05    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSigType ;
                *ReactionIndex = 9;
        }else if (  randReaction >= ReactionProb[0] && randReaction < ReactionProb[1]  ){
                //--reaction 10 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[1] = 0.10    Eth = 35 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = 0;
                *EmitParticle = iSiCl2gType ;
                *ReactionIndex = 10;
         }else if (  randReaction >= ReactionProb[1] && randReaction< ReactionProb[2]  ){
                //--reaction 11 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[2] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = 0;
                *EmitParticle = iSiCl2gType ;
                *ReactionIndex = 11;
        }else if (  randReaction >= ReactionProb[2] && randReaction < ReactionProb[3] ){
                //--reaction 12 : SiCl2(s) + Cl+ --> SiCl2(g) + Cl*        p0_ClIonReaction[3] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl2gType ;
                *ReactionIndex = 12;
        }else if(  randReaction >= ReactionProb[3] && randReaction < ReactionProb[4]  ){
                //--reaction 13 : SiCl3(s) + Cl+ --> SiCl3(g) + Cl*   p0_ClIonReaction[4] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl3gType ;
                *ReactionIndex = 13;
        }else{
                *ReflectedParticle = iClIonType;
                *EmitParticle = 0;
                *ReactionIndex = 0;
        }
        if ( cell::dNumMaterial[*iTag] == 0 )      cell::setStatus(*iTag, iVacuumStat, 1);
}



void cell::Cl2IonReaction(int* iTag, double* energy, double* IncidentAngle, double* ReactionProb, int* number_of_reactions,
                                                    int* ReflectedParticle, int* EmitParticle, int* ReactionIndex){

        random_device rd;
        default_random_engine generator( rd() );  //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randReaction = unif(generator);
        if (  randReaction < ReactionProb[0] ){
                //--reaction  15 : Si(s) + Cl2^+ --> Si(g) + Cl2*        p0_Cl2IonReaction[0] = 0.02    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSigType ;
                *ReactionIndex = 15;
        }else if ( randReaction >= ReactionProb[0] && randReaction < ReactionProb[1]  ){
                //--reaction 16 : SiCl(s) + Cl2^+ --> SiCl2(g) + Cl*        p0_Cl2IonReaction[1] = 0.20    Eth 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl2gType ;
                *ReactionIndex = 16;
        }else if (  randReaction >= ReactionProb[1] && randReaction < ReactionProb[2]  ){
                //--reaction 17 : SiCl2(s) + Cl2^+ --> SiCl2(g) + Cl2*        p0_Cl2IonReaction[2] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSiCl2gType ;
                *ReactionIndex = 17;
        }else if (  randReaction >= ReactionProb[2] && randReaction < ReactionProb[3]  ){
                //--reaction 18 : SiCl2(s) + Cl2^+ --> SiCl3(g) + Cl*        p0_Cl2IonReaction[3] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl3gType ;
                *ReactionIndex = 18;
         }else if (  randReaction >= ReactionProb[3] && randReaction < ReactionProb[4]  ){
                //--reaction 19 : SiCl3(s) + Cl2+ --> SiCl3(g) + Cl2*        p0_Cl2IonReaction[4] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSiCl3gType ;
                *ReactionIndex = 19;
        }else if (  randReaction >= ReactionProb[4] && randReaction < ReactionProb[5] ){
                //--reaction 20 : SiCl3(s) + Cl2+ --> SiCl4(g) + Cl*        p0_Cl2IonReaction[5] = 0.25    Eth =10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = 0 ;
                *ReactionIndex = 20;
        }else{
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = 0;
                *ReactionIndex = 0;
        }
        if ( cell::dNumMaterial[*iTag] == 0 )   cell::setStatus(*iTag, iVacuumStat, 1);
}


void cell::ArIonReaction(int* iTag, double* energy, double* IncidentAngle, double* ReactionProb, int* number_of_reactions,
                                                    int* ReflectedParticle, int* EmitParticle, int* ReactionIndex){
        random_device rd;
        std::default_random_engine generator( rd() );  //--random number generator
        std::uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        double randReaction = unif(generator);
        if (  randReaction < ReactionProb[0] && dNumSiClxs[*iTag][0] > 0 ){
                //--reaction 22 : Si(s) + Ar+ --> Si(g) + Ar*    p0_ArIonReaction[0] = 0.05    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iArIonType ;
                *EmitParticle = iSigType ;
                *ReactionIndex = 22;
        }else if (  randReaction >= ReactionProb[0] && randReaction < ReactionProb[1]  ){
                //--reaction 23 : SiCl(s) + Ar+ --> SiCl(g) + Ar*        p0_ArIonReaction[1] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiClgType;
                *ReactionIndex = 23;
        }else if (  randReaction >= ReactionProb[1] && randReaction < ReactionProb[2]  ){
                //--reaction 24 : SiCl2(s) + Ar+ --> SiCl2(g) + Ar*        p0_ArIonReaction[2] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiCl2gType;
                *ReactionIndex = 24;
        }else if (  randReaction >= ReactionProb[2] && randReaction < ReactionProb[3]  ){
                //--reaction 25 : SiCl3(s) + Ar+ --> SiCl3(g) + Ar*        p0_ArIonReaction[3] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[*iTag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[*iTag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiCl3gType;
                *ReactionIndex = 25;
        }else{
                *ReflectedParticle = iArIonType;
                *EmitParticle = 0;
                *ReactionIndex = 0;
        }
        if ( cell::dNumMaterial[*iTag] == 0 )   cell::setStatus(*iTag, iVacuumStat, 1);

}



//====================================================================================//


void cell::ElasticReflection(double* normReflectedVelocity, double* speed, double* Vel)
{
        Vel[X_dir] = *speed*normReflectedVelocity[X_dir];
        Vel[Y_dir] = *speed*normReflectedVelocity[Y_dir];
        Vel[Z_dir] = *speed*normReflectedVelocity[Z_dir];

}


void cell::InelasticReflection(double* normReflectedVelocity, double* GrazingAngle, double* mass, double* Vel, double* speed, double* energy)
{
        double EnergyDependentScalingFactor;
        double AngleDependentScalingFactor;
        double energy_in_unit_eV = *energy*Joule_to_eV;
        double epsilon_0 = 0.0;  //--unit: eV
        double epsilon_s = 50;  //--unit: eV
        double theta_0 = 30;  //--unit: degree
        double gamma_0 = 0.85;  //--a scaling factor, no unit
        if ( energy_in_unit_eV < epsilon_0 ){
                EnergyDependentScalingFactor = 0;
        }else if ( energy_in_unit_eV >= epsilon_0 && energy_in_unit_eV <= epsilon_s ){
                EnergyDependentScalingFactor = (energy_in_unit_eV - epsilon_0)/(epsilon_s - epsilon_0);
        }else if (  energy_in_unit_eV > epsilon_s){
                EnergyDependentScalingFactor = 1;
        }

        if( *GrazingAngle  < theta_0){
                AngleDependentScalingFactor = (theta_0 - *GrazingAngle)/(theta_0);
        }else if (*GrazingAngle >= theta_0){
                AngleDependentScalingFactor = 0.0;
        }
        //cout << "energy = " << *energy << endl;
        //cout << "speed = " << *speed<< endl;
        //cout << AngleDependentScalingFactor << endl;
        *energy = *energy * gamma_0 * AngleDependentScalingFactor * EnergyDependentScalingFactor; //--updated energy
        *speed = sqrt(2* *energy/(*mass));                                                                                                                                              //--updated speed
        //cout << "energy = " << *energy << endl;
        //cout << "speed = " << *speed<< endl;

        Vel[X_dir] = *speed*normReflectedVelocity[X_dir];                                                                                                     //--updated x component of velocity
        Vel[Y_dir] = *speed*normReflectedVelocity[Y_dir];                                                                                                     //--updated y component of velocity
        Vel[Z_dir] = *speed*normReflectedVelocity[Z_dir];                                                                                                     //--updated z component of
}



void cell::DiffusiveReflection(double* Temperature, double* normSurfaceNormal, double* ReemissionCosineLawPower,
                                                               double* mass, double* Vel, double* speed, double* energy)
{
        random_device rd;
        std::default_random_engine generator( rd() );
        std::normal_distribution<double> distribution(0.0, 1.0);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        double randTheta = unif(generator);
        double randPhi = unif(generator);
        double theta = acos(   pow( randTheta, 1/(*ReemissionCosineLawPower+1) )   );
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
                double v_mp = sqrt(2*BoltzmannConstant_SI*(*Temperature)/(*mass));  //--the most probable velocity = sqrt(2kT/m)
                Vel[X_dir] = v_mp*distribution(generator);
                Vel[Y_dir] = v_mp*distribution(generator);
                Vel[Z_dir] = v_mp*distribution(generator);
                *speed = sqrt(Vel[X_dir]*Vel[X_dir]+Vel[Y_dir]*Vel[Y_dir]+Vel[Z_dir]*Vel[Z_dir]);
        }while(*speed == 0);
        *energy = 0.5*(*mass)*(*speed)*(*speed);


        //--modify particle reflected velocity by new speed, theta, and phi
        //--Use the formula that I derived
        Vel[X_dir] = *speed*(   cos(alpha)*cos(beta)*sin(theta)*cos(phi)+cos(alpha)*sin(beta)*cos(theta)-sin(alpha)*sin(theta)*sin(phi)   );
        Vel[Y_dir] = *speed*(   sin(alpha)*cos(beta)*sin(theta)*cos(phi)+sin(alpha)*sin(beta)*cos(theta)+cos(alpha)*sin(theta)*sin(phi)   );
        Vel[Z_dir] = *speed*( -sin(beta)*sin(theta)*cos(phi)+cos(beta)*cos(theta)  );
}



