//
//  MainPrimaryGeneratorAction.c
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#include "MainPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

MainPrimaryGeneratorAction* MainPrimaryGeneratorAction::fgInstance = 0;

const MainPrimaryGeneratorAction* MainPrimaryGeneratorAction::Instance()
{
    // Static acces function via G4RunManager
    
    return fgInstance;
}

MainPrimaryGeneratorAction::MainPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
{
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    
    // default particle kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="proton");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1,0.,-1));
    //Need to fix this
    fParticleGun->SetParticleEnergy(6.*MeV);
    
    fgInstance = this;
}

MainPrimaryGeneratorAction::~MainPrimaryGeneratorAction()
{
    delete fParticleGun;
    fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MainPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //this function is called at the begining of ecah event
    //
    
    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get Envelope volume
    // from G4LogicalVolumeStore.
    
    G4double envSizeXY = 0;
    G4double envSizeZ = 0;
    G4LogicalVolume* envLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    G4Box* envBox = NULL;
    if ( envLV ) envBox = dynamic_cast<G4Box*>(envLV->GetSolid());
    if ( envBox ) {
        envSizeXY = envBox->GetXHalfLength()*2.;
        envSizeZ = envBox->GetZHalfLength()*2.;
    }
    else  {
        G4cerr << "Envelope volume of box shape not found." << G4endl;
        G4cerr << "Perhaps you have changed geometry." << G4endl;
        G4cerr << "The gun will be place in the center." << G4endl;
    }
    
    G4double size = 0.8;
    G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
    G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
    G4double z0 = -0.5 * envSizeZ;
    
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
