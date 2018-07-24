// Standard Includes:                                                                                                 
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

// Root Includes
#include "TApplication.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"

// Generator
#include "Generator.hh"

// Addition QoL libraries
#include "PhysicalConstants.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"

int main(int argc, char **argv)
{

  Generator *compton = new Generator();

  compton->SetBeamEnergy(5);                                   // Default beam energy
  compton->SetLaserEnergy(2.33e-9);                            // Laser energy (eV)
  compton->SetPolarization(compton->fPolarization);            // Defined as P=-P_o (left), P=P_o (right)

  compton->GetOptions(argv);
  if(argc == 0){
    Sys::SysError << __FUNCTION__ << " >>> Must include arguments! Use --help to see command list." << Sys::endl;
    exit(1);
  }
  if(!(compton->fHaloGenerator) && !(compton->fComptonGenerator)){ 
    Sys::SysError << __FUNCTION__ << " >>> Must set generator type flag." << Sys::endl;
    exit(1);
  }
  Sys::SysCout << "Polarization: " << compton->fPolarization << Sys::endl;
  if(compton->fGraphicsShow) compton->InitGraphicsEngine(argc, argv); 

  compton->Initialize();

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  compton->GenerateAsymmetry((char *)"");                     // The char * casting removes a deprecatred warning caused by difference between char * in C and C++
  compton->GetFunction((char *)"asym")->Draw();

  std::cout << 2*TMath::Pi()*compton->GetFunction((char *)"cs")->Integral(0,1) << std::endl;

  canvas->SaveAs("theory.C");

  gRandom->SetSeed(0);
  
  Sys::SysCout << "<<<< Random seed: " << gRandom->GetSeed() << Sys::endl;

  TF1 *beamx = new TF1("beamx", Generator::BeamEnvelope, -compton->upper_limit, compton->upper_limit, 1); // positions in cm the beam envelope is centered around zero
  beamx->SetParameter(0, compton->sigma_x);                  // positions in cm
  beamx->SetNpx(10000);

  TF1 *beamy = new TF1("beamy", Generator::BeamEnvelope, -compton->upper_limit, compton->upper_limit, 1); // positions in cm the beam envelope is centered around zero
  beamy->SetParameter(0, compton->sigma_y);                  // positions in cm
  beamy->SetNpx(10000);

  TF1 *halox;
  TF1 *haloy;

  if(compton->fHaloGenerator){

    halox = new TF1("halox", Generator::HaloFunctionX, compton->cutoffx, compton->upper_limit, 2); // positions in cm
    halox->SetParameters(compton->sigma_x, compton->halo_scale_x); // positions in cm
    halox->SetNpx(100000);

    haloy = new TF1("haloy", Generator::HaloFunctionY, compton->cutoffy, compton->upper_limit, 2); // positions in cm
    haloy->SetParameters(compton->sigma_y, compton->halo_scale_y); // positions in cm
    haloy->SetNpx(100000);
  }

  compton->OpenOutputFile();

  int sign = 1;

  for(int i = 0; i < (int)(compton->GetNumberEvents()); i++)
    {
      sign *= -1;    // Alternate sign of randomly sampled beam position. This is more efficient than sampling from a wider distribution.
      
      if(compton->fComptonGenerator) compton->SetEventVertex(-29.29464 + beamx->GetRandom(-7.0, 7.0), 
       							     beamy->GetRandom(-7.0, 7.0), 
       							     -2287.855); // add in the gaussian nature of the electron beam
      
      
      if(compton->fHaloGenerator) compton->SetEventVertex(-29.29464 + beamx->GetRandom(-7.0, 7.0) + sign*halox->GetRandom(-7.0, 7.0), 
							  beamy->GetRandom(-7.0, 7.0) + sign*haloy->GetRandom(-7.0, 7.0), 
							  -2187.855);    // add in the gaussian nature of the electron beam
      
      if(i % 10000 == 0) Sys::SysCout << i << Sys::endl;
      if(compton->fComptonGenerator){
	compton->CalculateKinematics();
	compton->ProcessComptonEvent();
      }
      if(compton->fHaloGenerator){
	compton->CalculateKinematics();
	compton->ProcessHaloEvent();
      }
    }
  
  compton->CloseOutputFile();
  
  if(compton->fGraphicsShow){
    compton->RunGraphicsEngine(); 
  }
  
  Sys::SysCout << "Finished." << Sys::endl;
  
  
  return 0;
}
