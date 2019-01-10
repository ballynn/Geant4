//
//  MainRunAction.c
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#include "MainRunAction.hh"
#include "MainPrimaryGeneratorAction.hh"
#include "MainEventAction.hh"
#include "MainSteppingAction.hh"
    // use of other actions
    // - primary generator: to get info for printing about the primary
    // - event action: to get and reset accumulated energy sums
    // - stepping action: to get info about accounting volume
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

MainRunAction::MainRunAction()
: G4UserRunAction()
{
    // add new units for dose
    //
    const G4double milligray = 1.e-3*gray;
    const G4double microgray = 1.e-6*gray;
    const G4double nanogray  = 1.e-9*gray;
    const G4double picogray  = 1.e-12*gray;
    
    new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
    new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
    new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
    new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);
}

MainRunAction::~MainRunAction()
{}

void MainRunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    
    //inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
    //initialize event cumulative quantities
    MainEventAction::Instance()->Reset();
}

void MainRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4int nofEvents = aRun->GetNumberOfEvent();
    if (nofEvents == 0) return;
    
    // Compute dose
    //
    G4double energySum  = MainEventAction::Instance()->GetEnergySum();
    G4double energy2Sum = MainEventAction::Instance()->GetEnergy2Sum();
    G4double rms = energy2Sum - energySum*energySum/nofEvents;
    if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;
    
    G4double mass = MainSteppingAction::Instance()->GetVolume()->GetMass();
    G4double dose = energySum/mass;
    G4double rmsDose = rms/mass;
    
    
// Run conditions
//
const G4ParticleGun* particleGun
    = MainPrimaryGeneratorAction::Instance()->GetParticleGun();
G4String particleName
    = particleGun->GetParticleDefinition()->GetParticleName();
G4double particleEnergy = particleGun->GetParticleEnergy();
    
// Print
//
G4cout
<< "\n--------------------End of Run------------------------------\n"
<< " The run consists of " << nofEvents << " "<< particleName << " of "
<<   G4BestUnit(particleEnergy,"Energy")
<< "\n Dose in scoring volume "
<< MainSteppingAction::Instance()->GetVolume()->GetName() << " : "
<< G4BestUnit(dose,"Dose")
<< " +- "                   << G4BestUnit(rmsDose,"Dose")
<< "\n------------------------------------------------------------\n"
<< G4endl;
}

