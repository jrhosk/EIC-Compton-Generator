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


int main(int argc, char **argv)
{

  TApplication theApp("App", &argc, argv);

  int nevents = 10;

  char *filename;

  filename = argv[1];

  Generator *compton = new Generator();

  compton->SetBeamEnergy(5);
  compton->SetLaserEnergy(2.33e-9);
  compton->SetPolarization(1.0);

  compton->Initialize();
  //  compton->GetFunction()->SetNpx(1000);

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();
  compton->GetFunction()->Draw("surf1");

  gRandom->SetSeed(0);
  
  std::cout << "<<< Random seed: " << gRandom->GetSeed() << std::endl;

  compton->SetEventVertex(-29.29464, 0.0, -2287.855);
  compton->OpenOutputFile(filename);

  for(int i = 0; i < nevents; i++)
    {
      compton->CalculateKinematics();
      compton->ProcessEvent();
    
    }
  compton->CloseOutputFile();

  std::cout << "Finished." << std::endl;

  theApp.Run();

  return 0;
}
