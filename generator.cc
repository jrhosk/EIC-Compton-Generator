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
#include "FontColor.hh"


int main(int argc, char **argv)
{

  Generator *compton = new Generator();

  compton->SetBeamEnergy(5);
  compton->SetLaserEnergy(2.33e-9);
  compton->SetPolarization(compton->fPolarization);            // Defined as P=-1 (left), P=1 (right)

  compton->GetOptions(argv);
  std::cout << "Polrization: " << compton->fPolarization << std::endl;
  if(compton->fGraphicsShow) compton->InitGraphicsEngine(argc, argv); 

  compton->Initialize();

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  compton->GenerateAsymmetry((char *)""); // The char * casting removes a deprecatred warning caused by difference between char * in C and C++
  compton->GetFunction((char *)"asym")->Draw();

  std::cout << 2*TMath::Pi()*compton->GetFunction((char *)"cs")->Integral(0,1) << std::endl;

  canvas->SaveAs("theory.C");

  gRandom->SetSeed(0);
  
  std::cout << green << "<<<< Random seed: " << gRandom->GetSeed() << white << std::endl;

  TF1 *beamx = new TF1("beamx", Generator::BeamEnvelope, -7.0, 7.0, 1); // positions in cm the beam envelope is centered around zero
  beamx->SetParameter(0, compton->sigma_x); // positions in cm
  beamx->SetNpx(10000);

  TF1 *beamy = new TF1("beamy", Generator::BeamEnvelope, -7.0, 7.0, 1); // positions in cm the beam envelope is centered around zero
  beamy->SetParameter(0, compton->sigma_y); // positions in cm
  beamy->SetNpx(10000);


  compton->OpenOutputFile();

  for(int i = 0; i < (int)(compton->GetNumberEvents()); i++)
    {

      compton->SetEventVertex(-29.29464 + beamx->GetRandom(-7.0, 7.0), 0.0 + beamy->GetRandom(-7.0, 7.0), -2287.855); // add in the gaussian nature of the electron beam

      if(i % 10000 == 0) std::cout << i << std::endl;
      compton->CalculateKinematics();
      compton->ProcessEvent();
    
    }

  compton->CloseOutputFile();

  if(compton->fGraphicsShow){
    compton->RunGraphicsEngine(); 
  }

  std::cout << green << "Finished." << white << std::endl;

  
  return 0;
}
