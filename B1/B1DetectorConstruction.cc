//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//LYNN EDIT
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 100*m, env_sizeZ = 20*m;
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
    
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  
    //
    // Envelope
    //
  G4Box*solidEnv =
    new G4Box("Envelope",
              env_sizeXY, env_sizeXY, env_sizeZ);
  G4LogicalVolume* logicenv =
    new G4LogicalVolume(solidEnv,
                        world_mat,
                        "Envelope");
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicenv,
                      "Envelope",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    
  //
  // Shape1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_WATER");
  G4ThreeVector pos1 = G4ThreeVector(0, 0, -30);
  G4double cyl_rad = 5*m;
  G4double cyl_height = 5*m;
  G4Tubs* solidcyl =
    new G4Tubs("Tissue",                       //its name
        0, cyl_rad, cyl_height, 0, 360*degree);  //its size
      
  G4LogicalVolume* Tissue =
    new G4LogicalVolume(solidcyl,            //its solid
                        shape1_mat,          //its material
                        "Tissue");         //its name
    
   new G4PVPlacement(0,                       //no rotation
                    pos1, //at (0,0,0)
                    Tissue,                //its logical volume
                    "Tissue",                //its name
                    logicenv,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    //
    // Bottom
    //
   G4Material* matbottom = nist->FindOrBuildMaterial("G4_Al");
   G4ThreeVector pos2 = G4ThreeVector(0, 0, 40);
   G4double xylength = 10*cm;
   G4double height = 5*cm;
   G4Box* bottom =
   new G4Box("Shape2",                       //its name
             xylength, xylength, height);//its size
    
    G4LogicalVolume* logicbottom =
    new G4LogicalVolume(bottom,            //its solid
                        matbottom,          //its material
                      "Bottom");         //its name
    
   new G4PVPlacement(0,                       //no rotation
                      pos2,                   //at (0,0,0)
                     logicbottom,                //its logical volume
                     "Bottom",              //its name
                     logicenv,              //its mother  volume
                     false,                   //no boolean operation
                     0,                       //copy number
                     checkOverlaps);          //overlaps checking
    
    
  // Hallbach shape
    G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_F");
    G4ThreeVector pos3 = G4ThreeVector(0, 0, -200);
    G4double tube_outrad = 10*cm;
    G4double tube_innerrad = 3*cm;
    G4double tube_height = 3*cm;

    //Magnetic Field
    G4UniformMagField* magField =
    new G4UniformMagField(G4ThreeVector(0,0,20.*tesla));
    
    G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldManager->SetDetectorField(magField);
    fieldManager->CreateChordFinder(magField);
    
    //Set Equation of Motion in Magnetic Field
    
    G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(magField);
    G4MagIntegratorStepper* pStepper = new G4ClassicalRK4(fEquation);
    G4MagInt_Driver*pIntgrDriver = new G4MagInt_Driver(1e-5*mm, pStepper, pStepper->GetNumberOfVariables());
    G4ChordFinder* pChordFinder = new G4ChordFinder(pIntgrDriver);
    pChordFinder->SetDeltaChord(1e-3*mm); //Maximum Miss Distance
    fieldManager->SetChordFinder(pChordFinder);
    fieldManager->SetDeltaOneStep(1e-3*mm);
    fieldManager->SetDeltaIntersection(1e-4*mm);
    
    //G4PropagatorInField* fieldPropogator = G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
     //fieldPropogator->SetMinimumEpsilonStep(1e-5*mm);
     //fieldPropogator->SetMaximuilomEpsilonStep(1e-2*mm);
     //fieldPropogator->SetLargestAcceptableStep(10*mm);
    
    G4Tubs* solidhalbach =
    new G4Tubs("Halbach",                       //its name
               tube_innerrad, tube_outrad, tube_height, 0, 360*degree);//its size
    
    G4LogicalVolume* logichalbach =
    new G4LogicalVolume(solidhalbach,            //its solid
                        shape2_mat,          //its material
                        "Halbach",          //its name
                        fieldManager);
    
    new G4PVPlacement(0,                       //no rotation
                      pos3,                   //at (0,0,0)
                      logichalbach,                //its logical volume
                      "Halbach",              //its name
                      logicenv,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    logichalbach->SetFieldManager(fieldManager, true);
    
    
  // Set Shape2 as scoring volume
  //
  fScoringVolume = Tissue;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
