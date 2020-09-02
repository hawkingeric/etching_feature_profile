#include "particle.h"
#include "etching.h"
#include "rand.h"

using namespace std;



//--Probability density function of Maxwell Boltzmann distribution
double Maxwell_Boltzmann_pdf(double mass, double T, double v){
        double v_p = sqrt(2*BoltzmannConstant_SI*T/mass); // the most probable speed
        return 4*PI*v*v/pow(PI, 1.5)/pow(v_p, 3.0)*exp( -pow(v/v_p, 2.0)  );
}


//--Cumulative density function of Maxwell Boltzmann distribution
double Maxwell_Boltzmann_cdf(double mass, double T, double v){
        double v_p = sqrt(2*BoltzmannConstant_SI*T/mass); // the most probable speed
        return erf(v/v_p) - 2/sqrt(PI)*(v/v_p)*exp( -pow(v/v_p, 2.0)  );
}


//--Calculate velocity of Cumulative density function of Maxwell Boltzmann distribution
double calc_v_cut(double mass, double T){
        double criterion = 0.0001;
        double dv = 0.1;
        double v_rms = sqrt(3*BoltzmannConstant_SI*T/mass);  // the root mean square velocity = sqrt(3kT/m)
        double v = v_rms;
        while(  (1 - Maxwell_Boltzmann_cdf(mass, T, v)  ) > criterion){
                v = v + dv;
        };
        return v;
}


void particle::setInitialPosition(double Lx, double Ly, double Lz, double dx, double dy, double dz)
{
        random_device rd;
        default_random_engine generator( rd() ); //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution
        particle::dPos[X_dir] = unif(generator)*Lx;
        particle::dPos[Y_dir] = unif(generator)*Ly;
        particle::dPos[Z_dir] = Lz-dz*0.5;
        particle::iPos[X_dir] = int(particle::dPos[X_dir]/dx);
        particle::iPos[Y_dir] = int(particle::dPos[Y_dir]/dy);
        particle::iPos[Z_dir] = int(particle::dPos[Z_dir]/dz);
}


void particle::setInitialType(vector<double>& ParticleProb_for_incident_particle){
        random_device rd;
        default_random_engine generator( rd() ); //--random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  //--random number uniform distribution

        int number_of_incident_particle = ParticleProb_for_incident_particle.size() ;
        double CumuProb [number_of_incident_particle] = {0.0};
        CumuProb[0] = ParticleProb_for_incident_particle[0];

        for(int i = 1; i < number_of_incident_particle; i ++){
                CumuProb[i] = CumuProb[i-1] + ParticleProb_for_incident_particle[i] ;
        }

        double RandParticle = unif(generator);
        if (RandParticle >= 0 && RandParticle < CumuProb[0]){
                particle::ParticleType = iClRadicalType;
        }else if ( RandParticle >= CumuProb[0] && RandParticle < CumuProb[1]){
                particle::ParticleType = iSigType;
        }else if ( RandParticle >= CumuProb[1] && RandParticle < CumuProb[2]){
                particle::ParticleType = iSiClgType;
        }else if ( RandParticle >= CumuProb[2] && RandParticle < CumuProb[3]){
                particle::ParticleType = iSiCl2gType;
        }else if ( RandParticle >= CumuProb[3] && RandParticle < CumuProb[4]){
                 particle::ParticleType = iSiCl3gType;
        }else if ( RandParticle >= CumuProb[4] && RandParticle < CumuProb[5]){
                 particle::ParticleType = iClIonType;
        }else if ( RandParticle >= CumuProb[5] && RandParticle < CumuProb[6]){
                 particle::ParticleType = iCl2IonType;
        }else if ( RandParticle >= CumuProb[6] && RandParticle < CumuProb[7] ){
                particle::ParticleType = iArIonType;
        }else{
                particle::ParticleType = 0;
        }

}


double particle::setInitialSpeed(double temperature, double mass, double v_cut)
{
        double v_mp = sqrt(2*BoltzmannConstant_SI*temperature/mass);  //--the most probable velocity = sqrt(2kT/m)
        double p_max = Maxwell_Boltzmann_pdf(mass, temperature, v_mp);  //--the maximal probability = p(v_mp)
        double v_rand, p_rand, prob_density;
        random_device rd;
        default_random_engine generator( rd() ); //--Random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  // --Random number uniform distribution
        do
        {
                v_rand = unif(generator)*v_cut;          //--Generate first random numbr/
                prob_density = Maxwell_Boltzmann_pdf(mass, temperature, v_rand);   //--Mawwell Boltzmann distribution probability density funcion
                p_rand = unif(generator)*p_max;         //--Generate second random number
        }while (prob_density < p_rand);
        return v_rand;
}


double particle::setInitialTheta(int type)
{
         double incident_angle;
         double center_angle = PI;
         random_device rd;
         default_random_engine generator( rd() ); // random number generator
         uniform_real_distribution<double> unif(0.0, 1.0);  // random number uniform distribution
         if(type == iArIonType || type == iClIonType || type == iCl2IonType ){  //--For ion particles
                 double sigma = 0.001;
                 double p_max = 1/sigma/sqrt(2*PI);
                 double t_rand, p_rand, prob_density;
                 double scale = 2.35482; //FWHM = Full width at half maximum = 2 * sqrt( 2 * ln(2)  )
                 do{
                         t_rand = unif(generator)*scale*sigma - scale*sigma + center_angle;
                         prob_density = 1/sigma/sqrt(2*PI)*exp(-0.5*pow(  (t_rand-center_angle)/sigma,2.0 )  ); //--Gaussian distribution probability density function
                         p_rand = unif(generator)*p_max;
                 }while (prob_density < p_rand);
                 incident_angle = t_rand;
         }else if (  type == iClRadicalType ){   //--For particles withour charge
                 incident_angle = unif(generator)*PI/2+PI/2;          //--generate the random number between PI/2 to PI for theta
         }
         return incident_angle;
}


double particle::setReemitTheta(int type)
{
        double reemit_angle;
        random_device rd;
        default_random_engine generator( rd() ); // random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  // random number uniform distribution

        double sigma = 0.5;
        double center_angle = 0;
        double p_max = 1/sigma/sqrt(2*PI);
        double t_rand, p_rand, prob_density;
        double scale = 2.35482; //FWHM = Full width at half maximum = 2 * sqrt( 2 * ln(2)  )
        do{
                t_rand = center_angle - scale*sigma + unif(generator)*scale*sigma ;
                prob_density = 1/sigma/sqrt(2*PI)*exp(-0.5*pow(  (t_rand-center_angle)/sigma, 2.0 )  ) ; //--Gaussian distribution probability density function
                p_rand = unif(generator)*p_max ;
        }while (prob_density < p_rand) ;
        reemit_angle = t_rand ;

        return reemit_angle;
}


void particle::reflected_velocity_with_new_energy(double* norm_reflected_V, double* grazing_angle, double* reflected_velocity){
        double energy_factor;
        double angle_factor;
        double epsilon_0 = 0.0/Joule_to_eV;
        double epsilon_s = 50/Joule_to_eV;
        double theta_0 = 30;
        double gamma_0 = 0.85;

        if ( particle::energy < epsilon_0 ){
                energy_factor = 0;
        }else if ( particle::energy >= epsilon_0 && particle::energy <= epsilon_s ){
                energy_factor = (particle::energy - epsilon_0)/(epsilon_s - epsilon_0);
        }else if (particle::energy > epsilon_s){
                energy_factor = 1;
        }

        if( *grazing_angle  < theta_0){
                angle_factor = (theta_0 - *grazing_angle)/theta_0;
        }else if (*grazing_angle >= theta_0){
                angle_factor = 0.0;
        }

        particle::energy = gamma_0 * angle_factor * energy_factor * particle::energy; //--updated energy
        particle::speed = sqrt(2*particle::energy/particle::mass); //--updated speed

        particle::Vel[X_dir] = particle::speed*norm_reflected_V[X_dir];
        particle::Vel[Y_dir] = particle::speed*norm_reflected_V[Y_dir];
        particle::Vel[Z_dir] = particle::speed*norm_reflected_V[Z_dir];
}


