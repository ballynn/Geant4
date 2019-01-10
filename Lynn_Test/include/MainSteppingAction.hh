//
//  MainSteppingAction.h
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#ifndef MainSteppingAction_h
#define MainSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;

/// Stepping action class
///
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the
/// SetVolume() function. The accumulated energy deposit is reset for each
/// new event via the Reset() function from the event action.

class MainSteppingAction : public G4UserSteppingAction
{
    public:
    MainSteppingAction();
    virtual ~MainSteppingAction();
    
    // static access method
    static MainSteppingAction* Instance();
    
    // method from the base class
    virtual void UserSteppingAction(const G4Step*);
    
    // reset accumulated energy
    void Reset();
    
    // set methods
    void SetVolume(G4LogicalVolume* volume) { fVolume = volume; }
    
    // get methods
    G4LogicalVolume* GetVolume() const { return fVolume; }
    G4double GetEnergy() const { return fEnergy; }
    
    private:
    static MainSteppingAction* fgInstance;
    
    G4LogicalVolume* fVolume;
    G4double  fEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
