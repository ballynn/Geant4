//
//  MainEventAction.h
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#ifndef MainEventAction_h
#define MainEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class MainSteppingAction;

/// Event action class
///
/// It holds data member fEnergySum and fEnergy2Sum for accumulating
/// the event energy deposit its square event by event.
/// These data are then used in the run action to compute the dose.
/// The accumulated energy and enrgy square sums are reset for each
/// new run via the Reset() function from the run action.

class MainEventAction : public G4UserEventAction
{
    public:
    MainEventAction();
    virtual ~MainEventAction();
    
    // static access method
    static MainEventAction* Instance();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    void Reset();
    
    // get methods
    G4double GetEnergySum() const { return fEnergySum; }
    G4double GetEnergy2Sum() const { return fEnergy2Sum; }
    
    private:
    static MainEventAction* fgInstance;
    
    G4int     fPrintModulo;
    G4double  fEnergySum;
    G4double  fEnergy2Sum;
};


#endif

