#define Generator_cxx
#include "../include/Generator.hh"

#include "TF2.h"
#include "TMath.h"

#ifdef Generator_cxx

Generator::Generator(char *options){}

// Generator::Generator(double beam_e, double laser_e, double polar, char *options)
// {

//   beam_energy = beam_e;
//   laser_energy = laser_e;
//   polarization = polar;

// }

Generator::~Generator(){}


Double_t Generator::CrossSection(Double_t *x = 0, Double_t *par = 0)
{

  Double_t rho        = x[0];
  Double_t phi        = x[1];
  Double_t b_energy   = par[0];
  Double_t l_energy   = par[1];
  Double_t P          = par[2];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  // Unpolarized cross section

  Double_t term1 = (rho*rho*(1-alpha)*(1-alpha))/(1-rho*(1-alpha));
  Double_t term2 = TMath::Power((1-rho*(1+alpha))/(1-rho*(1-alpha)), 2);
  Double_t fdSig_dRho_0 = 2*TMath::Pi()*TMath::Power(electron_radius, 2)*alpha*(term1 + 1 + term2);

  // Polarized longitudinal cross section

  Double_t term3 = (1-rho*(1+alpha))*(1-( 1/(1-rho*(1-alpha))));
  Double_t fdSig_dRho_1 =2*TMath::Pi()*TMath::Power(electron_radius, 2)*alpha*term3;

  // Total cross section for zero transverse polarization

  Double_t fdSig_dRho = fdSig_dRho_0 + P*fdSig_dRho_1;

  return fdSig_dRho;

}

void Generator::Initialize(){

  cs = new TF2("cs", this, &Generator::CrossSection, 0.0, 1.0, 0.0, TMath::Pi(), 3, "Generator", "CrossSection");
  cs->SetParameters(beam_energy, laser_energy, polarization);
  cs->SetNpx(1000);
}

void Generator::SetBeamEnergy(double energy){ 
  beam_energy = energy;
}
void Generator::SetLaserEnergy(double energy){ 
  laser_energy = energy;
}
void Generator::SetLaserWaveLength(double lambda){ 
  if(lambda > 0) laser_energy = h_planck*c_light/lambda;
}
void Generator::SetPolarization(double polar){ 
  polar = polarization;
}

double Generator::GetBeamEnergy(){return beam_energy;}

double Generator::GetLaserEnergy(){return laser_energy;}

double Generator::GetPolarization(){return polarization;}

double Generator::GetElectronEnergy(){return kinematics.electron_energy;}

double Generator::GetPhotonEnergy(){return kinematics.photon_energy;}

double Generator::GetElectronTheta(){return kinematics.electron_theta;}

double Generator::GetPhotonTheta(){return kinematics.photon_theta;}

double Generator::GetElectronPx(){return kinematics.px;}

double Generator::GetElectronPy(){return kinematics.py;}

double Generator::GetElectronPz(){return kinematics.pz;}

double Generator::GetPhotonPx(){return kinematics.kx;}

double Generator::GetPhotonPy(){return kinematics.ky;}

double Generator::GetPhotonPz(){return kinematics.kz;}


void Generator::CalculateKinematics()
{
  double rho  = 0;

  cs->GetRandom2(rho, kinematics.photon_phi);

  kinematics.kmax = 4*alpha*laser_energy*std::pow(beam_energy/electron_mass_c2, 2); // The maximum scattered photon energy or minimum electron energy

  kinematics.photon_energy = rho*kinematics.kmax;                               
  kinematics.electron_phi = -kinematics.photon_phi;                    

  kinematics.electron_energy = std::sqrt(std::pow(beam_energy, 2) - std::pow(electron_mass_c2, 2)); 
  kinematics.photon_theta = std::sqrt(4*laser_energy*alpha/kinematics.photon_energy-std::pow(electron_mass_c2/beam_energy, 2));
  kinematics.electron_theta = std::asin(kinematics.photon_energy*std::sin(kinematics.photon_theta)/kinematics.electron_energy); // check this

  kinematics.px = kinematics.electron_energy*std::sin(kinematics.electron_theta)*std::sin(kinematics.electron_phi);
  kinematics.py = kinematics.electron_energy*std::sin(kinematics.electron_theta)*std::cos(kinematics.electron_phi);
  kinematics.py = kinematics.electron_energy*std::cos(kinematics.electron_theta);

  kinematics.kx = kinematics.photon_energy*std::sin(kinematics.photon_theta)*std::sin(kinematics.photon_phi);
  kinematics.ky = kinematics.photon_energy*std::sin(kinematics.photon_theta)*std::cos(kinematics.photon_phi);
  kinematics.ky = kinematics.photon_energy*std::cos(kinematics.photon_theta);

}

void Generator::CalculateKinematics(event *kinematic){}

void Generator::OpenOutputFile(char *filename)
{
  output.open(filename, std::fstream::out);

  if(!(output.is_open())){
    std::cerr << "Failure to open output file. Exiting." << std::endl;
    exit(1);
  }

}

void Generator::WriteHeader()
{
  event_counter++;

  output << "2 "
	 << beam_energy << " "
	 << laser_energy
	 << " 0." << " "      // Integrated cross section                                                                                                            
	 << polarization << " "
	 << "0. 0. 0. 0. 0. \n";

}
void Generator::WriteEvent(int index, int pid, double px, double py, double pz, double energy)
{

  output << index
	 << " 0. 1 "
	 << pid
	 << " 0 0 "
	 << px << " "
	 << py << " "
	 << pz << " "
	 << energy << " "
	 << "0." << " "
	 << kinematics.vx << " "
	 << kinematics.vy << " "
	 << kinematics.vz << "\n";
  
}

void Generator::ProcessEvent()
{

  WriteHeader(); // Write event header for a given event

  // Write event info for each final state particle.
  WriteEvent(1, pid_electron, kinematics.px, kinematics.py, kinematics.pz, kinematics.electron_energy);
  WriteEvent(2, pid_photon, kinematics.kx, kinematics.ky, kinematics.kz, kinematics.photon_energy);


}

void Generator::SetEventVertex(double x, double y, double z)
{

  kinematics.vx = x;    // Set scattering vertex coordinates (cm)
  kinematics.vy = y;
  kinematics.vz = z;

}

void Generator::CloseOutputFile()
{
  event_counter = 0;
  output.close();
}

TF2* Generator::GetFunction()
{
  return cs;
}

#endif
