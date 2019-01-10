//
//  MainEventAction.c
//  
//
//  Created by Lynn Sargent on 1/10/19.
//
#include "MainEventAction.hh"

#include "MainRunAction.hh"
#include "MainSteppingAction.hh"
// use of stepping action to get and reset accumulated energy

#include "G4RunManager.hh"
#include "G4Event.hh"

MainEventAction* MainEventAction::fgInstance = 0;

MainEventAction* MainEventAction::Instance()
{
    // Static acces function via G4RunManager
    
    return fgInstance;
}



MainEventAction::MainEventAction()
: G4UserEventAction(),
fPrintModulo(100),
fEnergySum(0.),
fEnergy2Sum(0.)
{
    fgInstance = this;
}


MainEventAction::~MainEventAction()
{
    fgInstance = 0;
}


void MainEventAction::BeginOfEventAction(const G4Event* event)
{
    G4int eventNb = event->GetEventID();
    if (eventNb%fPrintModulo == 0) {
        G4cout << "\n---> Begin of event: " << eventNb << G4endl;
    }
    
    // Reset accounted energy in stepping action
    MainSteppingAction::Instance()->Reset();
}

void MainEventAction::EndOfEventAction(const G4Event* /*event*/)
{
    // accumulate statistics
    G4double energy = MainSteppingAction::Instance()->GetEnergy();
    fEnergySum  += energy;
    fEnergy2Sum += energy*energy;
}

void MainEventAction::Reset()
{
    //reset cumulative quantities
    //
    fEnergySum = 0.;
    fEnergy2Sum = 0.;
}
