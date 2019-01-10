//
//  MainSteppingAction.c
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#include "MainSteppingAction.hh"

#include "MainDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

MainSteppingAction* MainSteppingAction::fgInstance = 0;

MainSteppingAction* MainSteppingAction::Instance()
{
    // Static acces function via G4RunManager
    
    return fgInstance;
}

MainSteppingAction::MainSteppingAction()
: G4UserSteppingAction(),
fVolume(0),
fEnergy(0.)
{
    fgInstance = this;
}

MainSteppingAction::~MainSteppingAction()
{
    fgInstance = 0;
}

void MainSteppingAction::UserSteppingAction(const G4Step* step)
{
    // get volume of the current step
    G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
    
    // check if we are in scoring volume
    if (volume != fVolume ) return;
    
    // collect energy and track length step by step
    G4double edep = step->GetTotalEnergyDeposit();
    fEnergy += edep;
}

void MainSteppingAction::Reset()
{
    fEnergy = 0.;
}

