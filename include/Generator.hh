#include "TF2.h"
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

#ifndef Generator_h
#define Generator_h 

class Generator {

private:

  struct event {

    double electron_theta;   // Electron polar angle (outgoing)
    double electron_phi;     // Electron azimuthal angle (outgoing)
    double photon_theta;     // Photon polar angle   (outgoing)
    double photon_phi;       // Photon azimuthal angle   (outgoing)

    double electron_energy;  // Scattered electron energy (eV)
    double photon_energy;    // Scattered photon energy (eV)

    double kmax;             // Maximum scattered photon energy

    double px;               // Electron x-momentum 
    double py;               // Electron y-momentum 
    double pz;               // Electron z-momentum 
    
    double kx;               // Photon x-momentum
    double ky;               // Photon y-momentum
    double kz;               // Photon z-momentum

    double vx;               // Event vertex X
    double vy;               // Event vertex Y
    double vz;               // Event vertex Z

  };

  TF2 *cs;

  int event_counter;

public:

  static const double electron_mass_c2 = 0.5109989461;      // MeV/c^2
  static const double electron_radius  = 2.817e-13;         // Classical electron radius (cm)
  static const double h_planck = 6.626070040e-34;           // Planck constant (J*s)
  static const double c_light  = 2.99792458e8;              // Speed of light (m/s)

  static const int pid_electron = 11;
  static const int pid_photon   = 22;

  std::fstream output;

  event kinematics;        // Defined above and holds the scattering kinematics

  double alpha;            // Kinematic variable
  double laser_energy;     // Initial photon energy (eV)
  double beam_energy;      // Initial electron energy (eV)
  double polarization;     // Beam polarization

  Generator(char *options = 0);
  // Generator(double, double, double, char *options);
  ~Generator();

  Double_t CrossSection(Double_t *x, Double_t *par);

  void Initialize();
  void SetBeamEnergy(double);
  void SetLaserEnergy(double);
  void SetLaserWaveLength(double);
  void SetPolarization(double);

  double GetBeamEnergy();
  double GetLaserEnergy();
  double GetPolarization();

  double GetElectronEnergy();
  double GetPhotonEnergy();

  double GetElectronTheta();
  double GetPhotonTheta();

  double GetElectronPx();
  double GetElectronPy();
  double GetElectronPz();

  double GetPhotonPx();
  double GetPhotonPy();
  double GetPhotonPz();

  void SetEventVertex(double, double, double);
  void CalculateKinematics();
  void CalculateKinematics(event *kinematic);

  void OpenOutputFile(char *filename);

  void WriteHeader();
  void WriteEvent(int, int, double, double, double, double);
  void ProcessEvent();
  void CloseOutputFile();

  TF2* GetFunction();


};

#endif
