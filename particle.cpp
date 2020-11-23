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
double CalculateSpeedCutoff(double mass, double T){
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
                particle::ParticleType = iClIonType;
        }else if ( RandParticle >= CumuProb[1] && RandParticle < CumuProb[2]){
                 particle::ParticleType = iCl2IonType;
        }else if ( RandParticle >= CumuProb[2] && RandParticle < CumuProb[3]){
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


double particle::setInitialTheta_by_Gaussian(double sigma)
{
         random_device rd;
         default_random_engine generator( rd() );
         uniform_real_distribution<double> unif(0.0, 1.0);
         double incident_angle;
         double center_angle = PI;
         double p_max = 1/sigma/sqrt(2*PI);
         double t_rand, p_rand, prob_density;
         double scale = 2.35482; //--FWHM = Full width at half maximum = 2 * sqrt( 2 * ln(2)  )
         do{
                 t_rand = unif(generator)*scale*sigma - scale*sigma + center_angle;
                 //--Gaussian distribution probability density function
                 prob_density = 1/sigma/sqrt(2*PI)*exp(-0.5*pow(  (t_rand-center_angle)/sigma,2.0 )  );
                 p_rand = unif(generator)*p_max;
         }while (prob_density < p_rand);
         incident_angle = t_rand;
         return incident_angle;
}


double particle::setReemittedTheta(double n)
{
        double reemit_angle;
        random_device rd;
        default_random_engine generator( rd() ); // random number generator
        uniform_real_distribution<double> unif(0.0, 1.0);  // random number uniform distribution
        double t_rand = unif(generator);
        reemit_angle = acos(   pow( t_rand, 1/(n+1) )   );
        return reemit_angle;
}

//--input: normReflectedVelocity, GrazingAngle, epsilon_0, epsilon_s, theta_0, gamma_0, E, speed, velocity;  output: ReflectedVelocity
void particle::ReflectedWithNewEnergy(double* normReflectedVelocity, double* GrazingAngle, double* epsilon_0, double* epsilon_s,
                                                                                       double* theta_0, double* gamma_0){
        double EnergyDependentScalingFactor;
        double AngleDependentScalingFactor;
        double Energy = particle::energy*Joule_to_eV;
        double Mass = particle::mass;

        if ( Energy < *epsilon_0 ){
                EnergyDependentScalingFactor = 0;
        }else if ( Energy >= *epsilon_0 && Energy <= *epsilon_s ){
                EnergyDependentScalingFactor = (Energy - *epsilon_0)/(*epsilon_s - *epsilon_0);
        }else if (  Energy > *epsilon_s){
                EnergyDependentScalingFactor = 1;
        }

        if( *GrazingAngle  < *theta_0){
                AngleDependentScalingFactor = (*theta_0 - *GrazingAngle)/(*theta_0);
        }else if (*GrazingAngle >= *theta_0){
                AngleDependentScalingFactor = 0.0;
        }

        particle::energy = particle::energy*(*gamma_0)*AngleDependentScalingFactor*EnergyDependentScalingFactor; //--updated energy
        particle::speed = sqrt(2*particle::energy/particle::mass); //--updated speed
        particle::Vel[X_dir] = particle::speed*normReflectedVelocity[X_dir];  //--updated x component of velocity
        particle::Vel[Y_dir] = particle::speed*normReflectedVelocity[Y_dir];  //--updated y component of velocity
        particle::Vel[Z_dir] = particle::speed*normReflectedVelocity[Z_dir];  //--updated z component of velocity
}


void particle::ReemittedWithNewDirection(double* normSurfaceNormal, double speed, double theta, double phi){
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

        //--Use the formula that I derived
        particle::Vel[X_dir] = speed*(   cos(alpha)*cos(beta)*sin(theta)*cos(phi)+cos(alpha)*sin(beta)*cos(theta)-sin(alpha)*sin(theta)*sin(phi)   );
        particle::Vel[Y_dir] = speed*(   sin(alpha)*cos(beta)*sin(theta)*cos(phi)+sin(alpha)*sin(beta)*cos(theta)+cos(alpha)*sin(theta)*sin(phi)   );
        particle::Vel[Z_dir] = speed*( -sin(beta)*sin(theta)*cos(phi)+cos(beta)*cos(theta)  );

        /*
        if (normSurfaceNormal[Z_dir] < 0){
                cout <<  "normSurfaceNormal = " << normSurfaceNormal[X_dir] << " " << normSurfaceNormal[Y_dir] << " " << normSurfaceNormal[Z_dir]<<endl;
                cout << "position = " << particle::iPos[X_dir] << " " << particle::iPos[Y_dir] << " " << particle::iPos[Z_dir] << endl;
                cout << "type = " << particle::ParticleType << endl;
                cin >> alpha ;
        }
        */
        //cout << "norm = " << pow(normSurfaceNormal[X_dir],2) + pow(normSurfaceNormal[Y_dir],2) + pow(normSurfaceNormal[Z_dir],2)<<endl;
        //cin.get();



}
