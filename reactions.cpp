#include "etching.h"
#include "cell.h"

using namespace std;

void cell::ClRadicalReaction(vector<double>& p0_ClRadicalReaction, int itag, int iNumMaterial, int* ReactionExecution, int* reaction_index){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);
        int number_of_reactions = p0_ClRadicalReaction.size();

        //--calculation of total prob, and cumulative prob
        vector<double> ReactionProb;
        double sumReactionProb = 0;
        double denominator = cell::dNumSiClxs[itag][0] + cell::dNumSiClxs[itag][1] + cell::dNumSiClxs[itag][2] * 2 + cell::dNumSiClxs[itag][3] * 2;
        ReactionProb.push_back(   p0_ClRadicalReaction[0]*cell::dNumSiClxs[itag][0]/denominator   );
        ReactionProb.push_back(   p0_ClRadicalReaction[1]*cell::dNumSiClxs[itag][1]/denominator   );
        ReactionProb.push_back(   p0_ClRadicalReaction[2]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   p0_ClRadicalReaction[3]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   p0_ClRadicalReaction[4]*cell::dNumSiClxs[itag][3]/denominator   );
        ReactionProb.push_back(   p0_ClRadicalReaction[5]*cell::dNumSiClxs[itag][3]/denominator   );
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
                cell::dNumSiClxs[itag][0]--;
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]++;
                *ReactionExecution = 1;
                *reaction_index = 1;
        }else if  (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 2 : SiCl(s) + Cl --> SiCl2(s)        p0_ClRadicalReaction[1] = 0.40
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]--;
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]++;
                *ReactionExecution = 1;
                *reaction_index = 2;
        }else if  (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 3 : SiCl2(s) + Cl --> SiCl3(s)              p0_ClRadicalReaction[2] = 0.30
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]++;
                *ReactionExecution = 1;
                *reaction_index = 3;
        }else if ( happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]   ){
                //--reaction 4 : SiCl2(s) + Cl --> SiCl(s) + Cl2         p0_ClRadicalReaction[3] = 0.02
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]++;
                *ReactionExecution = 1;
                *reaction_index = 4;
        }else if  ( happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 5 : SiCl3(s) + Cl --> SiCl4(g)                p0_ClRadicalReaction[4] = 0.0001
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *reaction_index = 5;
        }else if ( happen_or_not >= ReactionProb[4] && happen_or_not < ReactionProb[5]  )  {
                //--reaction 6 : SiCl3(s) + Cl --> SiCl2(s) + Cl2    p0_ClRadicalReaction[5] = 0.08
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]++;
                *ReactionExecution = 1;
                *reaction_index = 6 ;
        }else{
                *ReactionExecution = 0;
                *reaction_index = 0;
        }
        if ( cell::dNumMaterial[itag] == 0 )   cell::setStatus(itag, iVacuumStat, 1);
}



void cell::redeposition(vector<double>& p0_redeposition, int AdsorbParticle, int itag, int* ReactionExecution, int* reaction_index){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);


        //--reaction 8 : M(s) + SiClx(g)  --> M(s) + SiClx(s)   p = 0.02   x=1,2,3
        if ( AdsorbParticle == iSiClgType ){
                if (   happen_or_not < p0_redeposition[0] ){
                        //--reaction 8-1 : M(s)  +   SiCl(g)  -->  M(s)  +  SiCl(s)        p0_redeposition[0] = 0.02
                        if ( cell::iStatus[itag] == iVacuumStat)   cell::iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[itag][1]++;
                        #pragma omp atomic
                        cell::dNumMaterial[itag]++;
                        *ReactionExecution = 1;
                        *reaction_index = 8;
                }
        }else if ( AdsorbParticle == iSiCl2gType ){
                if(  happen_or_not < p0_redeposition[1] ){
                        //--reaction 8-2 : M(s) +   SiCl2(g)  -->  M(s)  +  SiCl2(s)        p0_redeposition[1] = 0.02
                        if ( cell::iStatus[itag] == iVacuumStat  )   cell::iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[itag][2]++;
                        #pragma omp atomic
                        cell::dNumMaterial[itag]++;
                        *ReactionExecution = 1;
                        *reaction_index = 8;
                }
        }else if (  AdsorbParticle == iSiCl3gType ){
                if (   happen_or_not < p0_redeposition[2] ){
                        //--reaction 8-3 : M(s)  +  SiCl3(g)  -->  M(s)  +  SiCl3(s)        p0_redeposition[2] = 0.02
                        if ( cell::iStatus[itag] == iVacuumStat  )   cell::iStatus[itag] = iSubstrateStat;
                        #pragma omp atomic
                        cell::dNumSiClxs[itag][3]++;
                        #pragma omp atomic
                        cell::dNumMaterial[itag]++;
                        *ReactionExecution = 1;
                        *reaction_index = 8;
                }
        }
}




void cell::ClIonReaction(vector<double>& Eth_ClIonReaction, double* E0, vector<double>& p0_ClIonReaction, vector<int>& type_ClIonReaction,
                                                    vector<double>& phys_sputter_prob, vector<double>& chem_sputter_prob,
                                                    int itag, int iNumMaterial, double E, double incident_angle, int* ReactionExecution, int* EmitParticle, int* reaction_index){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);
        int number_of_reactions = p0_ClIonReaction.size();

        //--Calculation of prob of energy and prob of angle
        vector<double> prob_of_energy;
        vector<double> prob_of_angle;
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;
        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_ClIonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_ClIonReaction[i])/(*E0 - Eth_ClIonReaction[i]) );
                        prob_of_energy.push_back( p0_ClIonReaction[i]*EnhanceFactor );
                }else if (E < Eth_ClIonReaction[i]){
                        prob_of_energy.push_back(0.0);
                }

                if( type_ClIonReaction[i] == P_sputtering ){
                        angle_index = int(floor(incident_angle/phys_sputter_prob.size()));
                        angle_interval = 90.0/(phys_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle.push_back(   phys_sputter_prob[angle_index]    );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }else if (type_ClIonReaction[i] == C_sputtering){
                        angle_index = int(floor(incident_angle/chem_sputter_prob.size()));
                        angle_interval = 90.0/(chem_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle.push_back( chem_sputter_prob[angle_index]  );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }
        }

        /*
        if (  cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
               cell::iStatus[cell::iID_NBR[itag][22]] == iVacuumStat  ){
               for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
              for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
            for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }
        */


        //--calculation of total prob, and cumulative prob
        vector<double> ReactionProb;
        double sumReactionProb = 0;
        double denominator = cell::dNumSiClxs[itag][0] + cell::dNumSiClxs[itag][1]*2 + cell::dNumSiClxs[itag][2] + cell::dNumSiClxs[itag][3] ;
        ReactionProb.push_back(   prob_of_energy[0]*prob_of_angle[0]*cell::dNumSiClxs[itag][0]/denominator   );
        ReactionProb.push_back(   prob_of_energy[1]*prob_of_angle[1]*cell::dNumSiClxs[itag][1]/denominator   );
        ReactionProb.push_back(   prob_of_energy[2]*prob_of_angle[2]*cell::dNumSiClxs[itag][1]/denominator   );
        ReactionProb.push_back(   prob_of_energy[3]*prob_of_angle[3]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   prob_of_energy[4]*prob_of_angle[4]*cell::dNumSiClxs[itag][3]/denominator   );
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
                cell::dNumSiClxs[itag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSigType ;
                *reaction_index = 9;
        }else if (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 10 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[1] = 0.10    Eth = 35 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 10;
         }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 11 : SiCl(s) + Cl+ --> SiCl2(g)        p0_ClIonReaction[2] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution =1;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 11;
        }else if (  happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3] ){
                //--reaction 12 : SiCl2(s) + Cl+ --> SiCl2(g) + Cl*        p0_ClIonReaction[3] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 12;
        }else if(  happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 13 : SiCl3(s) + Cl+ --> SiCl3(g) + Cl*   p0_ClIonReaction[4] = 0.5    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 13;
        }else{
                *ReactionExecution = 0;
                *reaction_index = 0;
        }
        if ( cell::dNumMaterial[itag] == 0 )      cell::setStatus(itag, iVacuumStat, 1);

}



void cell::Cl2IonReaction(vector<double>& Eth_Cl2IonReaction, double* E0, vector<double>& p0_Cl2IonReaction, vector<int>& type_Cl2IonReaction,
                                                      vector<double>& phys_sputter_prob, vector<double>& chem_sputter_prob,
                                                       int itag, int iNumMaterial, double E,  double incident_angle, int* ReactionExecution, int* EmitParticle, int* reaction_index){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);
        int number_of_reactions = p0_Cl2IonReaction.size();


        //--Calculation of prob of energy and prob of angle
        vector<double> prob_of_energy;
        vector<double> prob_of_angle;
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;
        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_Cl2IonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_Cl2IonReaction[i])/(*E0 - Eth_Cl2IonReaction[i]) );
                        prob_of_energy.push_back( p0_Cl2IonReaction[i]*EnhanceFactor );
                }else if (E < Eth_Cl2IonReaction[i]){
                        prob_of_energy.push_back(0.0);
                }
                if( type_Cl2IonReaction[i] == P_sputtering ){
                        angle_index = int(floor(incident_angle/phys_sputter_prob.size()));
                        angle_interval = 90.0/(phys_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle.push_back(   phys_sputter_prob[angle_index]    );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }else if (type_Cl2IonReaction[i] == C_sputtering){
                        angle_index = int(floor(incident_angle/chem_sputter_prob.size()));
                        angle_interval = 90.0/(chem_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle.push_back( chem_sputter_prob[angle_index]  );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }
         }

        /*
        if (  cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
               cell::iStatus[cell::iID_NBR[itag][22]] == iVacuumStat  ){
               for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
              for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
            for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }
        */

        //--calculation of total prob, and cumulative prob
        vector<double> ReactionProb;
        double sumReactionProb = 0;
        double denominator = cell::dNumSiClxs[itag][0] + cell::dNumSiClxs[itag][1] + cell::dNumSiClxs[itag][2]*2 + cell::dNumSiClxs[itag][3]*2 ;
        ReactionProb.push_back(   prob_of_energy[0]*prob_of_angle[0]*cell::dNumSiClxs[itag][0]/denominator   );
        ReactionProb.push_back(   prob_of_energy[1]*prob_of_angle[1]*cell::dNumSiClxs[itag][1]/denominator   );
        ReactionProb.push_back(   prob_of_energy[2]*prob_of_angle[2]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   prob_of_energy[3]*prob_of_angle[3]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   prob_of_energy[4]*prob_of_angle[4]*cell::dNumSiClxs[itag][3]/denominator   );
        ReactionProb.push_back(   prob_of_energy[5]*prob_of_angle[5]*cell::dNumSiClxs[itag][3]/denominator   );
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
                cell::dNumSiClxs[itag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSigType ;
                *reaction_index = 15;
        }else if ( happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 16 : SiCl(s) + Cl2^+ --> SiCl2(g) + Cl*        p0_Cl2IonReaction[1] = 0.20    Eth 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 16;
        }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 17 : SiCl2(s) + Cl2^+ --> SiCl2(g) + Cl2*        p0_Cl2IonReaction[2] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl2gType ;
                *reaction_index = 17;
        }else if ( happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]  ){
                //--reaction 18 : SiCl2(s) + Cl2^+ --> SiCl3(g) + Cl*        p0_Cl2IonReaction[3] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 18;
         }else if (  happen_or_not >= ReactionProb[3] && happen_or_not < ReactionProb[4]  ){
                //--reaction 19 : SiCl3(s) + Cl2+ --> SiCl3(g) + Cl2*        p0_Cl2IonReaction[4] = 0.25    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl3gType ;
                *reaction_index = 19;
        }else if (  happen_or_not >= ReactionProb[4] && happen_or_not < ReactionProb[5] ){
                //--reaction 20 : SiCl3(s) + Cl2+ --> SiCl4(g) + Cl*        p0_Cl2IonReaction[5] = 0.25    Eth =10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl4gType ;
                *reaction_index = 20;
        }else{
                *ReactionExecution = 0;
                *reaction_index = 0;
        }
        if ( cell::dNumMaterial[itag] == 0 )   cell::setStatus(itag, iVacuumStat, 1);
}




void cell::ArIonReaction(vector<double>& Eth_ArIonReaction, double* E0, vector<double>& p0_ArIonReaction, vector<int>& type_ArIonReaction,
                                                    vector<double>& phys_sputter_prob, vector<double>& chem_sputter_prob,
                                                    int itag, int iNumMaterial, double E, double incident_angle, int* ReactionExecution, int* EmitParticle, int* reaction_index){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);
        int number_of_reactions = p0_ArIonReaction.size();

        //--Calculation of prob of energy and prob of angle
        vector<double> prob_of_energy;
        vector<double> prob_of_angle;
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        double EnhanceFactor;
        for (int i = 0; i < number_of_reactions; i++){
                if(E >= Eth_ArIonReaction[i]){
                        EnhanceFactor = sqrt( (E - Eth_ArIonReaction[i])/(*E0 - Eth_ArIonReaction[i]) );
                        prob_of_energy.push_back( p0_ArIonReaction[i]*EnhanceFactor );
                }else if (E < Eth_ArIonReaction[i]){
                        prob_of_energy.push_back(0.0);
                }
                if( type_ArIonReaction[i] == P_sputtering ){
                        angle_index = int(floor(incident_angle/phys_sputter_prob.size()));
                        angle_interval = 90.0/(phys_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if ( incident_angle == angle ){
                                prob_of_angle.push_back(   phys_sputter_prob[angle_index]    );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }else if (type_ArIonReaction[i] == C_sputtering){
                        angle_index = int(floor(incident_angle/chem_sputter_prob.size()));
                        angle_interval = 90.0/(chem_sputter_prob.size()-1);
                        angle = double(angle_index) * angle_interval;
                        next_angle = double(angle_index+1) * angle_interval;
                        if (  incident_angle == angle ){
                                prob_of_angle.push_back( chem_sputter_prob[angle_index]  );
                        }else{
                                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                                prob_of_angle.push_back( (chem_sputter_prob[angle_index+1]+chem_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio)  );
                        }
                }
         }

        /*
        if (  cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
               cell::iStatus[cell::iID_NBR[itag][22]] == iVacuumStat  ){
               for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
              for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }else if ( cell::iStatus[cell::iID_NBR[itag][12]] == iVacuumStat && cell::iStatus[cell::iID_NBR[itag][14]] == iVacuumStat &&
                          cell::iStatus[cell::iID_NBR[itag][4]] == iVacuumStat ){
            for (int i = 0; i < number_of_reactions; i++)  prob_of_angle[i] = 1;
        }
        */


        //--calculation of total prob, and cumulative prob
        vector<double> ReactionProb;
        double sumReactionProb = 0;
        double denominator = cell::dNumSiClxs[itag][0] + cell::dNumSiClxs[itag][1] + cell::dNumSiClxs[itag][2] + cell::dNumSiClxs[itag][3] ;
        ReactionProb.push_back(   prob_of_energy[0]*prob_of_angle[0]*cell::dNumSiClxs[itag][0]/denominator   );
        ReactionProb.push_back(   prob_of_energy[1]*prob_of_angle[1]*cell::dNumSiClxs[itag][1]/denominator   );
        ReactionProb.push_back(   prob_of_energy[2]*prob_of_angle[2]*cell::dNumSiClxs[itag][2]/denominator   );
        ReactionProb.push_back(   prob_of_energy[3]*prob_of_angle[3]*cell::dNumSiClxs[itag][3]/denominator   );
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


        if (  happen_or_not < ReactionProb[0] && cell::dNumSiClxs[itag][0] > 0 ){
                //--reaction 22 : Si(s) + Ar+ --> Si(g) + Ar*    p0_ArIonReaction[0] = 0.05    Eth = 25 eV    physical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][0]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSigType ;
                *reaction_index = 22;
        }else if (  happen_or_not >= ReactionProb[0] && happen_or_not < ReactionProb[1]  ){
                //--reaction 23 : SiCl(s) + Ar+ --> SiCl(g) + Ar*        p0_ArIonReaction[1] = 0.20    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][1]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiClgType;
                *reaction_index = 23;
        }else if (  happen_or_not >= ReactionProb[1] && happen_or_not < ReactionProb[2]  ){
                //--reaction 24 : SiCl2(s) + Ar+ --> SiCl2(g) + Ar*        p0_ArIonReaction[2] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][2]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl2gType;
                *reaction_index = 24;
        }else if (  happen_or_not >= ReactionProb[2] && happen_or_not < ReactionProb[3]  ){
                //--reaction 25 : SiCl3(s) + Ar+ --> SiCl3(g) + Ar*        p0_ArIonReaction[3] = 0.50    Eth = 10 eV    chemical sputtering
                #pragma omp atomic
                cell::dNumSiClxs[itag][3]--;
                #pragma omp atomic
                cell::dNumMaterial[itag]--;
                *ReactionExecution = 1;
                *EmitParticle = iSiCl3gType;
                *reaction_index = 25;
        }else{
                *ReactionExecution = 0;
                *reaction_index = 0;
        }
        if ( cell::dNumMaterial[itag] == 0 )   cell::setStatus(itag, iVacuumStat, 1);
}




void cell::IonMaskReaction(vector<double>& phys_sputter_prob, int itag, double E,  double incident_angle, int* ReactionExecution  ){
        random_device rd;
        default_random_engine generator( rd() );
        uniform_real_distribution<double> unif(0.0, 1.0);
        double happen_or_not = unif(generator);
        double prob_of_energy;
        double prob_of_angle;
        double EnhanceFactor;
        int angle_index;
        double angle_interval;
        double angle;
        double next_angle;
        double angle_ratio;
        if(E >= 15){
                EnhanceFactor = sqrt( (E - 15)/(100 - 15) );
                prob_of_energy = 0.01*EnhanceFactor ;
        }else if (E < 15){
                prob_of_energy = 0.0;
        }
        angle_index = int(floor(incident_angle/phys_sputter_prob.size()));
        angle_interval = 90.0/(phys_sputter_prob.size()-1);
        angle = double(angle_index) * angle_interval;
        next_angle = double(angle_index+1) * angle_interval;
        if ( incident_angle == angle ){
                prob_of_angle = phys_sputter_prob[angle_index];
        }else{
                angle_ratio = (next_angle - incident_angle)/(incident_angle - angle);
                prob_of_angle = (phys_sputter_prob[angle_index+1]+phys_sputter_prob[angle_index]*angle_ratio)/(1+angle_ratio) ;
        }


        if (happen_or_not <= prob_of_angle*prob_of_energy){
                #pragma omp atomic
                cell::dNumMask[itag]--;
                *ReactionExecution = 1;
        }else{
                *ReactionExecution = 0;
        }
        if (cell::dNumMask[itag] == 0)          cell::setStatus(itag, iVacuumStat, 1);
}



