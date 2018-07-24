#define Generator_cxx
#include "../include/Generator.hh"
#include "TF2.h"
#include "TF1.h"
#include "TMath.h"

// Custom Libraries
#include "PhysicalConstants.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"

#ifdef Generator_cxx

Generator::Generator(char *options)
{
  fFileName = (char *)"generator_lund.dat";
  fPolarization = 0.97;
  sigma_x = 226.6e-4;
  sigma_y = 99e-4;
  cutoffx = 0.0;
  cutoffy = 0.0;
  halo_scale_x = 10;
  halo_scale_y = 10;
  fGraphicsShow = false;
  fNumberEvents = 1000;
  fHaloGenerator = false;
  fComptonGenerator = false;
  upper_limit = 7.3024;
}

Generator::~Generator(){}


double Generator::CrossSection(double *x = 0, double *par = 0)
{

  double rho        = x[0];
  double b_energy   = par[0];
  double l_energy   = par[1];
  double P          = par[2];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  // Unpolarized cross section

  double term1 = (rho*rho*(1-alpha)*(1-alpha))/(1-rho*(1-alpha));
  double term2 = TMath::Power((1-rho*(1+alpha))/(1-rho*(1-alpha)), 2);
  double fdSig_dRho_0 = TMath::Power(electron_radius, 2)*alpha*(term1 + 1 + term2);

  // Polarized longitudinal cross section

  double term3 = (1-rho*(1+alpha))*(1-(1/TMath::Power( (1-rho*(1-alpha)),2) ));
  double fdSig_dRho_1 =TMath::Power(electron_radius, 2)*alpha*term3;

  // Total cross section for zero transverse polarization

  double fdSig_dRho = fdSig_dRho_0 + P*fdSig_dRho_1;

  return fdSig_dRho;

}

double Generator::BeamEnvelope(double *x, double *par){
  // Simple gaussian describing the core beam envelope.

  double pos = x[0];

  double sig = par[0];
  double diff_cross_section = 0;

  double top  = TMath::Power(pos, 2);
  double bot  = TMath::Power(sig, 2);

  diff_cross_section = TMath::Exp(-(top)/(2*bot) );


  return diff_cross_section;
}

void Generator::Initialize(){

  cs = new TF1("cs", this, &Generator::CrossSection, 0.0, 1.0, 3, "Generator", "CrossSection");
  cs->SetParameters(beam_energy, laser_energy, fPolarization);
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
  fPolarization = polar;
}

double Generator::GetBeamEnergy(){return beam_energy;}

double Generator::GetLaserEnergy(){return laser_energy;}

double Generator::GetPolarization(){return fPolarization;}

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
  double electron_prime = 0;

  kinematics.rho = cs->GetRandom();

  kinematics.photon_phi = gRandom->Uniform(0, 2*TMath::Pi()); // Sample from a uniform distribution of phi

  kinematics.kmax = 4*alpha*laser_energy*std::pow(beam_energy/electron_mass_c2, 2); // The maximum scattered photon energy or minimum electron energy

  kinematics.photon_momentum = kinematics.rho*kinematics.kmax;                               
  kinematics.electron_phi = -kinematics.photon_phi;                    

  electron_prime = beam_energy + laser_energy - kinematics.photon_momentum;

  kinematics.electron_momentum = std::sqrt(std::pow(electron_prime, 2) - std::pow(electron_mass_c2, 2)); 

  kinematics.photon_theta = std::sqrt(4*laser_energy*alpha/kinematics.photon_momentum-std::pow(electron_mass_c2/beam_energy, 2));
  kinematics.electron_theta = std::asin(kinematics.photon_momentum*std::sin(kinematics.photon_theta)/kinematics.electron_momentum); 

  kinematics.px = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::sin(kinematics.electron_phi);
  kinematics.py = kinematics.electron_momentum*std::sin(kinematics.electron_theta)*std::cos(kinematics.electron_phi);
  kinematics.pz = kinematics.electron_momentum*std::cos(kinematics.electron_theta);

  kinematics.kx = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::sin(kinematics.photon_phi);
  kinematics.ky = kinematics.photon_momentum*std::sin(kinematics.photon_theta)*std::cos(kinematics.photon_phi);
  kinematics.kz = kinematics.photon_momentum*std::cos(kinematics.photon_theta);

  kinematics.asymmetry = RhoToAsymmetry(beam_energy, laser_energy, kinematics.rho);
  // PrintAsymmetryInfo();
}

void Generator::PrintAsymmetryInfo()
{
  Sys::SysMsg << kinematics.rho << "\t"
	      << kinematics.asymmetry
	      << Sys::endl;
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
  double b_energy   = par[0];
  double l_energy   = par[1];

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  double minus = rho*(1-alpha);
  double plus  = rho*(1+alpha);
  double term1 = 1/( (std::pow(minus, 2)/(1-minus)) + 1 + std::pow((1-plus)/(1-minus),2));
  double term2 = 1-plus;
  double term3 = 1-(1/std::pow(1-minus,2));

  double asymmetry = fPolarization*term1*term2*term3;

  return asymmetry;
}

double Generator::RhoToAsymmetry(double b_energy = 0, double l_energy = 0, double rho = 0)
{

  alpha = 1/(1 + (4*l_energy*b_energy)/(electron_mass_c2*electron_mass_c2));  

  double minus = rho*(1-alpha);
  double plus  = rho*(1+alpha);
  double term1 = 1/( (std::pow(minus, 2)/(1-minus)) + 1 + std::pow((1-plus)/(1-minus),2));
  double term2 = 1-plus;
  double term3 = 1-(1/std::pow(1-minus,2));

  double asymmetry = fPolarization*term1*term2*term3;

  return asymmetry;
}

void Generator::OpenOutputFile()
{

  output.open(fFileName, std::fstream::out);

  if(!(output.is_open())){
    Sys::SysError << __FUNCTION__ << " Failure to open output file. Exiting." << Sys::endl;
    exit(1);
  }

}

void Generator::OpenOutputFile(char *filename)
{

  output.open(filename, std::fstream::out);

  if(!(output.is_open())){
    Sys::SysError << __FUNCTION__ << " Failure to open output file. Exiting." << Sys::endl;
    exit(1);
  }

}

void Generator::WriteHeader()
{
  event_counter++;

  output << fNumberofParticles << " "
         << beam_energy << " "
         << laser_energy << " "
         << kinematics.kmax << " "
	 << fPolarization << " "
	 << kinematics.rho << " "
	 << kinematics.asymmetry << " "
	 << "0. 0. 0. \n";


}
void Generator::WriteEvent(int index, int pid, double px, double py, double pz, double momentum)
{

  output << index
         << " 0. 1 "
         << pid << " "
	 << " 0 0 "
         << px << " "
         << py << " "
         << pz << " "
         << momentum << " "
	 << alpha << " "         // alpha << "-29.29464 0.0 -2287.855\n";                                                                                                                       
         << kinematics.vx << " "
         << kinematics.vy << " "
         << kinematics.vz << "\n";
  
}

void Generator::ProcessComptonEvent()
{

  WriteHeader(); // Write event header for a given event

  // Write event info for each final state particle.
  WriteEvent(1, pid_photon, -kinematics.kx, -kinematics.ky, -kinematics.kz, kinematics.photon_momentum);
  WriteEvent(2, pid_electron, -kinematics.px, -kinematics.py, -kinematics.pz, kinematics.electron_momentum);
}

void Generator::ProcessHaloEvent()
{

  WriteHeader(); // Write event header for a given event

  kinematics.electron_momentum = beam_energy;

  // Write event info for each final state particle.
  WriteEvent(1, pid_electron, 0.0, 0.0, -kinematics.electron_momentum, kinematics.electron_momentum);
}

void Generator::BuildGeneratedAsymmetryPlot()
{

}

void Generator::PrintEvent()
{

  Sys::SysCout << __FUNCTION__ << " \n=====================================\n" << Sys::endl;

  Sys::SysCout << __FUNCTION__ << " >>> \tElectron \n" << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron theta: " << kinematics.electron_theta << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron phi: " << kinematics.electron_phi << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron Momentum: " << kinematics.electron_momentum << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> ---------------------------     " << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron px: " << kinematics.px << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron py: " << kinematics.py << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Electron pz: " << kinematics.pz << Sys::endl;

  Sys::SysCout << __FUNCTION__ << " \n\n\n" << Sys::endl;

  Sys::SysCout << __FUNCTION__ << " >>> \tPhoton \n" << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon theta: " << kinematics.photon_theta << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon phi: " << kinematics.photon_phi << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon Momentum: " << kinematics.photon_momentum << Sys::endl;
  Sys::SysCout << __FUNCTION__ << "     ---------------------------     " << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon px: " << kinematics.kx << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon py: " << kinematics.ky << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon pz: " << kinematics.kz << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> Photon max: " << kinematics.kmax << Sys::endl;

  Sys::SysCout << __FUNCTION__ << " >>> =====================================\n" << Sys::endl;

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

int Generator::GetNumberEvents(){return fNumberEvents;}

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

void Generator::GetOptions(char **options)
{

  Int_t i = 0;

  TString flag;

  while(options[i] != NULL){
    flag = options[i];

    if(flag.CompareTo("--help", TString::kExact) == 0){
      flag.Clear();
      PrintHelp();
    }

    if(flag.CompareTo("--graphics", TString::kExact) == 0){
      flag.Clear();
      fGraphicsShow = true;

      Sys::SysMsg << __FUNCTION__  << "<<<< Initializing TApplication for plots.\t" << Sys::endl;
    }
    if(flag.CompareTo("--halo", TString::kExact) == 0){
      flag.Clear();
      fHaloGenerator = true;
      fNumberofParticles = 1;

      Sys::SysMsg << __FUNCTION__ << "<<<< Initializing halo generator.\t" << Sys::endl;
    }
    if(flag.CompareTo("--compton", TString::kExact) == 0){
      flag.Clear();
      fComptonGenerator = true;
      fNumberofParticles = 2;

      Sys::SysMsg << __FUNCTION__ << "<<<< Initializing compton generator.\t" << Sys::endl;
    }
    if(flag.CompareTo("--filename", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fFileName = options[i + 1];
      Sys::SysMsg << __FUNCTION__ << "<<<< Output file set to: " << fFileName << Sys::endl;
    }
    if(flag.CompareTo("--polarization", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fPolarization = atof(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Polarization set to: " << fPolarization << Sys::endl;
    }
    if(flag.CompareTo("--energy", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      beam_energy = atof(options[i + 1]);
    }
    if(flag.CompareTo("--sigmax", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      sigma_x = atof(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Beam X-dispersion set to: " << sigma_x << Sys::endl;
    }
    if(flag.CompareTo("--sigmay", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      sigma_y = atof(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Beam Y-dispersion set to: " << sigma_y << Sys::endl;
    }
    if(flag.CompareTo("--halo-scale-x", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      halo_scale_x = atof(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Halo X multiplier set to: " << halo_scale_x << Sys::endl;
    }
    if(flag.CompareTo("--halo-scale-y", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      halo_scale_y = atof(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Halo Y multiplier set to: " << halo_scale_y << Sys::endl;
    }
    if(flag.CompareTo("--cutoffx", TString::kExact) == 0){
      flag.Clear();
      std::cout << "<<<< Setting cutoff X:\t" << options[i+1] << std::endl;
      cutoffx = atof(options[i+1]);  
    }
    if(flag.CompareTo("--cutoffy", TString::kExact) == 0){
      flag.Clear();
      std::cout << "<<<< Setting cutoff Y:\t" << options[i+1] << std::endl;
      cutoffy = atof(options[i+1]);  
    }
    if(flag.CompareTo("--upper_limit", TString::kExact) == 0){
      flag.Clear();
      std::cout << "<<<< Setting upper limit:\t" << options[i+1] << std::endl;
      upper_limit = atof(options[i+1]);  
    }
    if(flag.CompareTo("--events", TString::kExact) == 0){
      std::string option(options[i+1]);
      flag.Clear();
      fNumberEvents = atoi(options[i + 1]);
      Sys::SysMsg << __FUNCTION__ << "<<<< Number of events generated: " << fNumberEvents << Sys::endl;
    }
    i++;
  }
  if(fComptonGenerator){
    if(beam_energy == 3){
      sigma_x = 136e-4;
      sigma_y = 56e-4;
    }
    if(beam_energy == 5){
      sigma_x = 226.6e-4;
      sigma_y = 99e-4;
    }
    if(beam_energy == 10){
      sigma_x = 434e-4;
      sigma_y = 199e-4;
    }
    Sys::SysMsg << __FUNCTION__ << "<<<< Beam energy set to: " << beam_energy 
	      << "\n\tsigma x: " << sigma_x 
	      << "\n\tsigma y: " << sigma_y 
	      << Sys::endl;
  }
  if(fHaloGenerator){
    if(beam_energy == 3){
      sigma_x = 136e-4;
      sigma_y = 56e-4;
    }
    if(beam_energy == 5){
      sigma_x = 356.6e-4;
      sigma_y = 115e-4;
    }
    if(beam_energy == 10){
      sigma_x = 709e-4;
      sigma_y = 229e-4;
    }
    Sys::SysMsg << __FUNCTION__ << "<<<< Beam energy set to: " << beam_energy 
	      << "\n\tsigma x: " << sigma_x 
	      << "\n\tsigma y: " << sigma_y 
	      << Sys::endl;
  }
}

void Generator::PrintHelp()
{

  Sys::SysCout << __FUNCTION__ << " >>> --halo                  [Required] Flag that that turns on halo generation." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --compton               [Required] Flag that that turns on compton generation." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --filename <name>       [Required] Set output filename." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --events <int>          [Required] Set number of events." << Sys::endl;

  Sys::SysCout << __FUNCTION__ << " >>> --polarization <double>     [Optional] Set polarization of generated events. Defaults to 97%." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --graphics                  [Optional] Flag that turns on graphical output." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --energy <double>           [Optional] Set beam energy. Defaults to 5 GeV." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --sigmax, -sigmay <double>  [Optional] Set core beam width in cm." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --halo-scale-x, --halo-scale-y <double> [Optional] Set multiplier that defines halo width compared to beam width. Defaults to x10." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --cutoffx, --cutoffy <double> [Optional] Set lower bound for events pulled from halo distribution(cm). Defaults to 0." << Sys::endl;
  Sys::SysCout << __FUNCTION__ << " >>> --upper_limit <double>        [Optional] Set beam energy. Defaults to 5 GeV.Set upper limit for events pulled from halo distribution (cm). Defaults to 7.3 cm" << Sys::endl;

}

double Generator::HaloFunctionX(double *x, double *par){

  const double A = 7.2e-5;

  double pos = x[0];

  double sig = par[0];
  double multiplier = par[1];
  double diff_cross_section = 0;
  double top  = TMath::Power(pos, 2);
  double bot_h = TMath::Power(multiplier*sig, 2);
 
  diff_cross_section = A*TMath::Exp(-(top)/(2*bot_h) );

  return(diff_cross_section);
}

double Generator::HaloFunctionY(double *x, double *par){

  const double A = 7.2e-5;

  double pos = x[0];

  double sig = par[0];
  double multiplier = par[1];
  double diff_cross_section = 0;
  double top  = TMath::Power(pos, 2);
  double bot_h = TMath::Power(multiplier*sig, 2);
 
  diff_cross_section = A*TMath::Exp(-(top)/(2*bot_h) );

  return(diff_cross_section);
}

void Generator::InitGraphicsEngine(int Argc, char **Argv)
{
  Sys::SysMsg << __FUNCTION__ << "<<<< Initialize Graphics Engine." << Sys::endl;
  app = new TApplication("App", &Argc, Argv);

}

void Generator::RunGraphicsEngine()
{
  Sys::SysMsg << __FUNCTION__ << "<<<< Running Graphics Engine." << Sys::endl;
  app->Run();
}

#endif
