//
//  MainPrimaryGeneratorAction.h
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#ifndef MainPrimaryGeneratorAction_h
#define MainPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class MainDetectorConstruction;

/// The primary generator action class with particle gum.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class MainPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
    MainPrimaryGeneratorAction();
    virtual ~MainPrimaryGeneratorAction();
    
    // static access method
    static const MainPrimaryGeneratorAction* Instance();
    
    // method from the base class
    virtual void GeneratePrimaries(G4Event*);
    
    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    
    private:
    static MainPrimaryGeneratorAction* fgInstance;
    
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

