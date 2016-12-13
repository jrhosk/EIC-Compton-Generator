// Standard Includes:                                                                                                 
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

// Root Includes
#include "TApplication.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom.h"

// Generator
#include "Generator.hh"
#include "FontColor.hh"


int main(int argc, char **argv)
{


  Generator *compton = new Generator();

  compton->GetOptions(argv);
  if(compton->fGraphicsShow) compton->InitGraphicsEngine(argc, argv); 

  compton->SetBeamEnergy(5);
  compton->SetLaserEnergy(2.33e-9);
  compton->SetPolarization(-1.0);

  compton->Initialize();

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  // compton->GenerateAsymmetry((char *)""); // The char * casting removes a deprecatred warning caused by difference between char * in C and C++
  // compton->GetFunction((char *)"asym")->Draw();
  compton->GetFunction((char *)"cs")->Draw();

  gRandom->SetSeed(0);
  
  std::cout << green << "<<<< Random seed: " << gRandom->GetSeed() << white << std::endl;

  compton->SetEventVertex(-29.29464, 0.0, -2287.855);
  // compton->SetEventVertex(0.00474491, -0.0061297, 8000);
  compton->OpenOutputFile();

  for(int i = 0; i < (int)(compton->GetNumberEvents()); i++)
    {
      if(i % 1000 == 0) std::cout << i << std::endl;
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
