
#define Generator_cxx
#include "../include/Generator.hh"

#include "TF2.h"
#include "TF1.h"
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
  // Double_t phi        = x[1];
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

  cs = new TF1("cs", this, &Generator::CrossSection, 0.0, 1.0, 3, "Generator", "CrossSection");
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
  polarization = polar;
}

double Generator::GetBeamEnergy(){return beam_energy;}

double Generator::GetLaserEnergy(){return laser_energy;}

double Generator::GetPolarization(){return polarization;}

double Generator::GetElectronMomentum(){return kinematics.electron_momentum;}

double Generator::GetPhotonMomentum(){return kinematics.photon_momentum;}

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

  rho = cs->GetRandom();

  std::cout << ">>>> Rho " << rho << std::endl;

  kinematics.kmax = 4*alpha*laser_energy*std::pow(beam_energy/electron_mass_c2, 2); // The maximum scattered photon energy or minimum electron energy

  kinematics.photon_momentum = rho*kinematics.kmax;                               
  kinematics.electron_phi = -kinematics.photon_phi;                    

  kinematics.electron_momentum = std::sqrt(std::pow(beam_energy, 2) - std::pow(electron_mass_c2, 2)); 
  // kinematics.photon_theta = std::sqrt(4*laser_energy*alpha/kinematics.photon_momentum-std::pow(electron_mass_c2/beam_energy, 2));
  kinematics.photon_theta = std::sqrt( 4.*kinematics.photon_momentum/kinematics.kmax - 1./(alpha*std::pow(beam_energy/electron_mass_c2, 2)));
  kinematics.electron_theta = std::asin(kinematics.photon_momentum*std::sin(kinematics.photon_theta)/kinematics.electron_momentum); // check this

  kinematics.px = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::sin(kinematics.electron_phi);
  kinematics.py = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::cos(kinematics.electron_phi);
  kinematics.pz = kinematics.electron_momentum*std::cos(kinematics.electron_theta);

  kinematics.kx = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::sin(kinematics.photon_phi);
  kinematics.ky = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::cos(kinematics.photon_phi);
  kinematics.kz = kinematics.photon_momentum*std::cos(kinematics.photon_theta);

}

void Generator::CalculateKinematics(event *kinematic){}

void Generator::GenerateAsymmetry(char* options)
{
  asym = new TF1("asym", this, &Generator::CalculateAsymmetry, 0.0, 1.0, 2, "Generator", "CalculateAsymmetry");
  asym->SetParameters(beam_energy, laser_energy);
  asym->SetNpx(1000);
}

double Generator::CalculateAsymmetry(double *x = 0, double *par = 0)
{

  double rho = x[0];
  Double_t b_energy   = par[0];
  Double_t l_energy   = par[1];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  double minus = rho*(1-alpha);
  double plus  = rho*(1+alpha);
  double term1 = 1/( (std::pow(minus, 2)/(1-minus)) + 1 + std::pow((1-plus)/(1-minus),2));
  double term2 = 1-plus;
  double term3 = 1-(1/std::pow(1-minus,2));

  double asymmetry = term1*term2*term3;

  return asymmetry;
}

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
         << kinematics.kmax << " "
	 << polarization << " "
         << "0. 0. 0. 0. 0. \n";


}
void Generator::WriteEvent(int index, int pid, double px, double py, double pz, double momentum)
{
  // output << index
  // 	 << " 0. 1 "
  // 	 << pid
  // 	 << " 0. 0. "
  // 	 << px << " "
  // 	 << py << " "
  // 	 << pz << " "
  // 	 << momentum << " "
  // 	 << "0." << " "
  // 	 << kinematics.vx << " "
  // 	 << kinematics.vy << " "
  // 	 << kinematics.vz << "\n";

  output << index
         << " 0. 1 "
         << pid
         << " 0 0 "
         << px << " "
         << py << " "
         << pz << " "
         << momentum << " "
         << alpha << "-29.29464 0.0 -2287.855\n";
  
}

void Generator::ProcessEvent()
{

  WriteHeader(); // Write event header for a given event

  // Write event info for each final state particle.
  WriteEvent(1, pid_photon, -kinematics.kx, -kinematics.ky, -kinematics.kz, kinematics.photon_momentum);
  WriteEvent(2, pid_electron, -kinematics.px, -kinematics.py, -kinematics.pz, kinematics.electron_momentum);

  //  PrintEvent();

}

void Generator::PrintEvent()
{

  std::cout << "\n=====================================\n" << std::endl;

  std::cout << "<<<< Electron \n" << std::endl;
  std::cout << "Electron theta: " << kinematics.electron_theta << std::endl;
  std::cout << "Electron phi: " << kinematics.electron_phi << std::endl;
  std::cout << "Electron Momentum: " << kinematics.electron_momentum << std::endl;
  std::cout << "     ---------------------------     " << std::endl;
  std::cout << "Electron px: " << kinematics.px << std::endl;
  std::cout << "Electron py: " << kinematics.py << std::endl;
  std::cout << "Electron pz: " << kinematics.pz << std::endl;

  std::cout << "\n\n\n" << std::endl;

  std::cout << "<<<< Photon \n" << std::endl;
  std::cout << "Photon theta: " << kinematics.photon_theta << std::endl;
  std::cout << "Photon phi: " << kinematics.photon_phi << std::endl;
  std::cout << "Photon Momentum: " << kinematics.photon_momentum << std::endl;
  std::cout << "     ---------------------------     " << std::endl;
  std::cout << "Photon px: " << kinematics.kx << std::endl;
  std::cout << "Photon py: " << kinematics.ky << std::endl;
  std::cout << "Photon pz: " << kinematics.kz << std::endl;
  std::cout << "\nPhoton max: " << kinematics.kmax << std::endl;

  std::cout << "\n=====================================\n" << std::endl;

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

TF1* Generator::GetFunction(char *option)
{

  std::string opt = std::string(option);

  if(strcmp(option, "cs") == 0)
    {
      return cs;
    }
    if(strcmp(option, "asym") == 0)
    {
      return asym;
    }
    return NULL;
}

#endif
