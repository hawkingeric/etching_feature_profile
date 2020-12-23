#include "etching.h"
#include "cell.h"
#include "rand.h"

using namespace std;

void cell::ClRadicalReaction(int itag, int* ReflectedParticle, int* reaction_index){

        double happen_or_not = ran3(&iRandTag);
        double p0_ClRadicalReaction [6] = {0.99, 0.40, 0.30, 0.02, 0.0001, 0.08};
        int number_of_reactions = 6;//p0_ClRadicalReaction.size();

        //--calculation of total prob, and cumulative prob
        double ReactionProb  [number_of_reactions];
        double sumReactionProb = 0;
        double denominator = dNumSiClxs[itag][0] + dNumSiClxs[itag][1] + dNumSiClxs[itag][2] * 2 + dNumSiClxs[itag][3] * 2;
        ReactionProb[0] = p0_ClRadicalReaction[0]*dNumSiClxs[itag][0]/denominator;
        ReactionProb[1] = p0_ClRadicalReaction[1]*dNumSiClxs[itag][1]/denominator;
        ReactionProb[2] = p0_ClRadicalReaction[2]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[3] = p0_ClRadicalReaction[3]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[4] = p0_ClRadicalReaction[4]*dNumSiClxs[itag][3]/denominator;
        ReactionProb[5] = p0_ClRadicalReaction[5]*dNumSiClxs[itag][3]/denominator;
        for( int i = 0; i < number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }
        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }


        if  ( happen_or_not < ReactionProb[0]  ){
                //--reaction 1 : Si(s) + Cl --> SiCl(s)        p0_ClRadicalReaction[0] = 0.99
                #pragma omp atomic
                dNumSiClxs[itag][0]--;
                #pragma omp atomic
                dNumSiClxs[itag][1]++;
                *ReflectedParticle = 0;
                *reaction_index = 1;
        }else if  (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 2 : SiCl(s) + Cl --> SiCl2(s)        p0_ClRadicalReaction[1] = 0.40
                #pragma omp atomic
                dNumSiClxs[itag][1]--;
                #pragma omp atomic
                dNumSiClxs[itag][2]++;
                *ReflectedParticle = 0;
                *reaction_index = 2;
        }else if  (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 3 : SiCl2(s) + Cl --> SiCl3(s)              p0_ClRadicalReaction[2] = 0.30
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumSiClxs[itag][3]++;
                *ReflectedParticle = 0;
                *reaction_index = 3;
        }else if ( happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]   ){
                //--reaction 4 : SiCl2(s) + Cl --> SiCl(s) + Cl2         p0_ClRadicalReaction[3] = 0.02
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumSiClxs[itag][1]++;
                *ReflectedParticle = 0;
                *reaction_index = 4;
        }else if  ( happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 5 : SiCl3(s) + Cl --> SiCl4(g)                p0_ClRadicalReaction[4] = 0.0001
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = 0;
                *reaction_index = 5;
        }else if ( happen_or_not >= ReactionProb[4] && happen_or_not < ReactionProb[5]  )  {
                //--reaction 6 : SiCl3(s) + Cl --> SiCl2(s) + Cl2    p0_ClRadicalReaction[5] = 0.08
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumSiClxs[itag][2]++;
                *ReflectedParticle = 0;
                *reaction_index = 6 ;
        }else{
                *ReflectedParticle = iClRadicalType;
                *reaction_index = 0;
        }
        if ( dNumMaterial[itag] == 0 ){
                iStatus[itag] = iVacuumStat;
        }
}


void cell::Redeposition( int AdsorbParticle, int itag, int* ReflectedParticle, int* reaction_index){

        double happen_or_not = ran3(&iRandTag);
        double p0_redeposition [3] = {0.02, 0.02, 0.02};

        //--reaction 8 : M(s) + SiClx(g)  --> M(s) + SiClx(s)   p = 0.02   x=1,2,3
        if ( AdsorbParticle == iSiClgType ){
                if (   happen_or_not < p0_redeposition[0] ){
                        //--reaction 8-1 : M(s)  +   SiCl(g)  -->  M(s)  +  SiCl(s)        p0_redeposition[0] = 0.02
                        if ( iStatus[itag] == iVacuumStat)   iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        dNumSiClxs[itag][1]++;
                        #pragma omp atomic
                        dNumMaterial[itag]++;
                        *ReflectedParticle = 0;
                        *reaction_index = 8;
                }else{
                        *ReflectedParticle = iSiClgType;
                }
        }else if ( AdsorbParticle == iSiCl2gType ){
                if(  happen_or_not < p0_redeposition[1] ){
                        //--reaction 8-2 : M(s) +   SiCl2(g)  -->  M(s)  +  SiCl2(s)        p0_redeposition[1] = 0.02
                        if ( iStatus[itag] == iVacuumStat  )   iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        dNumSiClxs[itag][2]++;
                        #pragma omp atomic
                        dNumMaterial[itag]++;
                        *ReflectedParticle = 0;
                        *reaction_index = 8;
                }else{
                        *ReflectedParticle = iSiCl2gType;
                }
        }else if (  AdsorbParticle == iSiCl3gType ){
                if (   happen_or_not < p0_redeposition[2] ){
                        //--reaction 8-3 : M(s)  +  SiCl3(g)  -->  M(s)  +  SiCl3(s)        p0_redeposition[2] = 0.02
                        if ( iStatus[itag] == iVacuumStat  )   iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        dNumSiClxs[itag][3]++;
                        #pragma omp atomic
                        dNumMaterial[itag]++;
                        *ReflectedParticle = 0;
                        *reaction_index = 8;
                }else{
                        *ReflectedParticle = iSiCl3gType;
                }
        }
}



void cell::ClIonReaction(int itag, double E, double incident_angle, int* ReflectedParticle, int* EmitParticle, int* reaction_index){

        double happen_or_not = ran3(&iRandTag);
        //cout << "happen_or_not in Cl ion reaction = " << happen_or_not << endl;
        double Eth_ClIonReaction [5] = {25, 35, 10, 10, 10};
        double E0 = 100;
        double p0_ClIonReaction [5] = {0.05, 0.10, 0.20, 0.50, 0.50};
        int type_ClIonReaction [5] = {0, 0, 1, 1, 1};
        double phys_sputter_prob [10] = {0.55, 0.555, 0.60, 0.555, 0.85, 1.3, 1.355, 1.0, 0.75, 0.0};
        double chem_sputter_prob [10] = {1.00, 1.000, 1.00, 1.000, 1.00, 0.9, 0.800, 0.6, 0.30, 0.0};
        int number_of_reactions = 5;//p0_ClIonReaction.size();

        //--Calculation of prob of energy and prob of angle
        double prob_of_energy [number_of_reactions];
        double prob_of_angle [number_of_reactions];
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;

        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_ClIonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_ClIonReaction[i])/(E0 - Eth_ClIonReaction[i]) );
                        prob_of_energy[i] = p0_ClIonReaction[i]*EnhanceFactor;
                }else if (E < Eth_ClIonReaction[i]){
                        prob_of_energy[i] = 0.0;
                }

                if( type_ClIonReaction[i] == P_sputtering ){
                        angle_index = int( floor(incident_angle/10)  );
                        angle_interval = 10.0;
                        angle = double(angle_index)*angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle[i] = phys_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }else if (type_ClIonReaction[i] == C_sputtering){
                        angle_index = int(floor(incident_angle/10));
                        angle_interval = 10.0;
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle[i] = chem_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }
        }

        //--calculation of total prob, and cumulative prob
        double ReactionProb [number_of_reactions];
        double sumReactionProb = 0;
        double denominator = dNumSiClxs[itag][0] + dNumSiClxs[itag][1]*2 + dNumSiClxs[itag][2] + dNumSiClxs[itag][3] ;
        ReactionProb[0] = prob_of_energy[0]*prob_of_angle[0]*dNumSiClxs[itag][0]/denominator;
        ReactionProb[1] = prob_of_energy[1]*prob_of_angle[1]*dNumSiClxs[itag][1]/denominator;
        ReactionProb[2] = prob_of_energy[2]*prob_of_angle[2]*dNumSiClxs[itag][1]/denominator;
        ReactionProb[3] = prob_of_energy[3]*prob_of_angle[3]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[4] = prob_of_energy[4]*prob_of_angle[4]*dNumSiClxs[itag][3]/denominator;
        for( int i = 0; i < number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }
        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }

        if ( happen_or_not < ReactionProb[0]   ){
                //--reaction 9 : Si(s) + Cl+ --> Si(g) + Cl*        p0_ClIonReaction[0] = 0.05    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][0]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSigType ;
                *reaction_index = 9;
        }else if (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 10 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[1] = 0.10    Eth = 35 eV    physical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][1]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = 0;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 10;
         }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 11 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[2] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][1]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = 0;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 11;
        }else if (  happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3] ){
                //--reaction 12 : SiCl2(s) + Cl+ --> SiCl2(g) + Cl*        p0_ClIonReaction[3] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 12;
        }else if(  happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 13 : SiCl3(s) + Cl+ --> SiCl3(g) + Cl*   p0_ClIonReaction[4] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 13;
        }else{
                *ReflectedParticle = iClIonType;
                *EmitParticle = 0;
                *reaction_index = 0;
        }
        if ( dNumMaterial[itag] == 0 ){
                //cout << "here" << endl;
                iStatus[itag] = iVacuumStat;
                dNumMaterial[itag] = 0;
        }
}



void cell::Cl2IonReaction(int itag, double E, double incident_angle, int* ReflectedParticle, int* EmitParticle, int* reaction_index){

        double happen_or_not = ran3(&iRandTag);
        //cout << "happen_or_not in Cl2 ion reaction = " << happen_or_not << endl;
        double Eth_Cl2IonReaction [6] = {25, 10, 10, 10, 10, 10};
        double E0 = 100;
        double p0_Cl2IonReaction [6] = {0.02, 0.20, 0.25, 0.25, 0.25, 0.25};
        int type_Cl2IonReaction [6] = {0, 1, 1, 1, 1, 1};
        int number_of_reactions = 6;//p0_Cl2IonReaction.size();
        double phys_sputter_prob [10] = {0.55, 0.555, 0.60, 0.555, 0.85, 1.3, 1.355, 1.0, 0.75, 0.0};
        double chem_sputter_prob [10] = {1.00, 1.000, 1.00, 1.000, 1.00, 0.9, 0.800, 0.6, 0.30, 0.0};

        //--Calculation of prob of energy and prob of angle
        double prob_of_energy [number_of_reactions];
        double prob_of_angle [number_of_reactions];
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;

        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_Cl2IonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_Cl2IonReaction[i])/(E0 - Eth_Cl2IonReaction[i]) );
                        prob_of_energy[i] = p0_Cl2IonReaction[i]*EnhanceFactor;
                }else if (E < Eth_Cl2IonReaction[i]){
                        prob_of_energy[i] = 0.0;
                }
                if( type_Cl2IonReaction[i] == P_sputtering ){
                        angle_index = int(  floor(incident_angle/10)  );
                        angle_interval = 10.0;
                        angle = double(angle_index)*angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle[i] = phys_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }else if (type_Cl2IonReaction[i] == C_sputtering){
                        angle_index = int(  floor(incident_angle/10)   );
                        angle_interval = 10.0;
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle[i] = chem_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }
         }

        //--calculation of total prob, and cumulative prob
        double ReactionProb [number_of_reactions];
        double sumReactionProb = 0;
        double denominator = dNumSiClxs[itag][0] + dNumSiClxs[itag][1] + dNumSiClxs[itag][2]*2 + dNumSiClxs[itag][3]*2 ;
        ReactionProb[0] = prob_of_energy[0]*prob_of_angle[0]*dNumSiClxs[itag][0]/denominator;
        ReactionProb[1] = prob_of_energy[1]*prob_of_angle[1]*dNumSiClxs[itag][1]/denominator;
        ReactionProb[2] = prob_of_energy[2]*prob_of_angle[2]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[3] = prob_of_energy[3]*prob_of_angle[3]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[4] = prob_of_energy[4]*prob_of_angle[4]*dNumSiClxs[itag][3]/denominator;
        ReactionProb[5] = prob_of_energy[5]*prob_of_angle[5]*dNumSiClxs[itag][3]/denominator;
        for( int i = 0; i < number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }
        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }


        if (  happen_or_not < ReactionProb[0] ){
                //--reaction  15 : Si(s) + Cl2^+ --> Si(g) + Cl2*        p0_Cl2IonReaction[0] = 0.02    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][0]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSigType ;
                *reaction_index = 15;
        }else if ( happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 16 : SiCl(s) + Cl2^+ --> SiCl2(g) + Cl*        p0_Cl2IonReaction[1] = 0.20    Eth 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][1]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 16;
        }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 17 : SiCl2(s) + Cl2^+ --> SiCl2(g) + Cl2*        p0_Cl2IonReaction[2] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 17;
        }else if ( happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]  ){
                //--reaction 18 : SiCl2(s) + Cl2^+ --> SiCl3(g) + Cl*        p0_Cl2IonReaction[3] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 18;
         }else if (  happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 19 : SiCl3(s) + Cl2+ --> SiCl3(g) + Cl2*        p0_Cl2IonReaction[4] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 19;
        }else if (  happen_or_not >= ReactionProb[4] && happen_or_not < ReactionProb[5] ){
                //--reaction 20 : SiCl3(s) + Cl2+ --> SiCl4(g) + Cl*        p0_Cl2IonReaction[5] = 0.25    Eth =10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iClIonType;
                *EmitParticle = 0 ;
                *reaction_index = 20;
        }else{
                *ReflectedParticle = iCl2IonType;
                *EmitParticle = 0;
                *reaction_index = 0;
        }
        if ( dNumMaterial[itag] == 0 ){
                iStatus[itag] = iVacuumStat;
                dNumMaterial[itag] = 0;
        }
}


void cell::ArIonReaction(int itag,  double E, double incident_angle, int* ReflectedParticle, int* EmitParticle, int* reaction_index){

        double happen_or_not = ran3(&iRandTag);
        //cout << "happen_or_not in Ar ion reaction = " << happen_or_not << endl;
        double Eth_ArIonReaction [4] = {25, 10, 10, 10};
        double E0 = 100;
        double p0_ArIonReaction [4] = {0.05, 0.20, 0.50, 0.50};
        int type_ArIonReaction [4] = {0, 1, 1, 1};
        int number_of_reactions = 4;//p0_ArIonReaction.size();
        double phys_sputter_prob [10] = {0.55, 0.555, 0.60, 0.555, 0.85, 1.3, 1.355, 1.0, 0.75, 0.0};
        double chem_sputter_prob [10] = {1.00, 1.000, 1.00, 1.000, 1.00, 0.9, 0.800, 0.6, 0.30, 0.0};

        //--Calculation of prob of energy and prob of angle
        double prob_of_energy [number_of_reactions];
        double prob_of_angle [number_of_reactions];
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;

        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_ArIonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_ArIonReaction[i])/(E0 - Eth_ArIonReaction[i]) );
                        prob_of_energy[i] = p0_ArIonReaction[i]*EnhanceFactor;
                }else if (E < Eth_ArIonReaction[i]){
                        prob_of_energy[i] = 0.0;
                }
                if( type_ArIonReaction[i] == P_sputtering ){
                        angle_index = int(  floor(incident_angle/10)   );
                        angle_interval = 10.0;
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle[i] = phys_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }else if (type_ArIonReaction[i] == C_sputtering){
                        angle_index = int(  floor(incident_angle/10)   );
                        angle_interval = 10.0;
                        angle = double(angle_index)*angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle[i] = chem_sputter_prob[angle_index];
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle[i] = (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio);
                        }
                }
         }


        //--calculation of total prob, and cumulative prob
        double ReactionProb [number_of_reactions];
        double sumReactionProb = 0;
        double denominator = dNumSiClxs[itag][0] + dNumSiClxs[itag][1] + dNumSiClxs[itag][2] + dNumSiClxs[itag][3] ;
        ReactionProb[0] = prob_of_energy[0]*prob_of_angle[0]*dNumSiClxs[itag][0]/denominator;
        ReactionProb[1] = prob_of_energy[1]*prob_of_angle[1]*dNumSiClxs[itag][1]/denominator;
        ReactionProb[2] = prob_of_energy[2]*prob_of_angle[2]*dNumSiClxs[itag][2]/denominator;
        ReactionProb[3] = prob_of_energy[3]*prob_of_angle[3]*dNumSiClxs[itag][3]/denominator;
        for( int i = 0; i < number_of_reactions; i++){
                sumReactionProb += ReactionProb[i];
        }
        if( sumReactionProb > 1.0){
                ReactionProb[0] /= sumReactionProb ;
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i]/sumReactionProb ;
                }
        }else{
                for( int i = 1; i < number_of_reactions; i++){
                        ReactionProb[i] = ReactionProb[i-1] + ReactionProb[i] ;
                }
        }


        if (  happen_or_not < ReactionProb[0] && dNumSiClxs[itag][0] > 0 ){
                //--reaction 22 : Si(s) + Ar+ --> Si(g) + Ar*    p0_ArIonReaction[0] = 0.05    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][0]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iArIonType ;
                *EmitParticle = iSigType ;
                *reaction_index = 22;
        }else if (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 23 : SiCl(s) + Ar+ --> SiCl(g) + Ar*        p0_ArIonReaction[1] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][1]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiClgType;
                *reaction_index = 23;
        }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 24 : SiCl2(s) + Ar+ --> SiCl2(g) + Ar*        p0_ArIonReaction[2] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][2]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiCl2gType;
                *reaction_index = 24;
        }else if (  happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]  ){
                //--reaction 25 : SiCl3(s) + Ar+ --> SiCl3(g) + Ar*        p0_ArIonReaction[3] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                dNumSiClxs[itag][3]--;
                #pragma omp atomic
                dNumMaterial[itag]--;
                *ReflectedParticle = iArIonType;
                *EmitParticle = iSiCl3gType;
                *reaction_index = 25;
        }else{
                *ReflectedParticle = iArIonType;
                *EmitParticle = 0;
                *reaction_index = 0;
        }
        if ( dNumMaterial[itag] == 0 ){
                iStatus[itag] = iVacuumStat;
                dNumMaterial[itag] = 0;
        }

}
