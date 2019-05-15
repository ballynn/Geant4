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
// $Id$
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
   // use of stepping action to set the accounting volume

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

#include "G4UniformGravityField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include "G4TransportationManager.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4RepleteEofM.hh"
#include "G4PropagatorInField.hh"
#include "G4UniformMagField.hh"


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
  G4double env_sizeXY = 20*m, env_sizeZ = 100*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
      //World Material
    
    G4double z, a, fractionmass, density, pressure, temperature;
    G4String name, symbol;
    G4int ncomponents;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
    a = 14.01*g/mole;
    G4Element*elN = new G4Element(name = "Nitrogen", symbol = "N", z = 7, a);
    
    a = 44.00*g/mole;
    G4Element*elCO2 = new G4Element(name = "Carbon Dioxide", symbol = "CO2", z = 8, a);
    
    a = 39.95*g/mole;
    G4Element*elAr = new G4Element(name = "Argon", symbol = "Ar", z = 9, a);
    
    a = 16.00*g/mole;
    G4Element*elO = new G4Element(name = "Oxygen", symbol = "O", z = 10, a);
    
    a = 28.01*g/mole;
    G4Element*elCO = new G4Element(name = "Carbon Monoxide", symbol  = "O", z = 11, a);
    
    density = 0.00002*g/cm3;   //0.020kg/m^3 from NASA mars fact sheet
    temperature = 200.15*kelvin;
    pressure = 0.6*atmosphere;
    //pressure = 0.00678*atmosphere;
    //pressure = 840*kg/(m*s);                                                        //reference A. Firan Paper
    G4Material* world_mat = new G4Material("Mars_Air", density, ncomponents = 5, kStateGas, temperature, pressure);
    world_mat->AddElement(elN, fractionmass = 1.89*perCent);
    world_mat->AddElement(elCO2, fractionmass = 96.0*perCent);
    world_mat->AddElement(elAr, fractionmass = 1.93*perCent);
    world_mat->AddElement(elO, fractionmass = 0.145*perCent);
    world_mat->AddElement(elCO, fractionmass = 0.1*perCent);
  
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
  
//  G4UniformGravityField* gravfield =
//    new G4UniformGravityField(G4ThreeVector(0.0, 0.0, -3.72076*CLHEP::m/CLHEP::s/CLHEP::s/CLHEP::c_light));
//
//  G4RepleteEofM* equation = new G4RepleteEofM(gravfield);
//
//  G4FieldManager* fieldManager =
//  G4TransportationManager::GetTransportationManager()->GetFieldManager();
//  fieldManager->SetDetectorField(gravfield);
//
//
////Equation of Motion
//    G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,8);
//    G4double minStep           = 0.01*mm;
//
//   G4ChordFinder* chordFinder =
//    new G4ChordFinder((G4MagneticField*)gravfield,minStep,stepper);
//
//    // Set accuracy parameters
//    G4double deltaChord        = 3.0*mm;
//    chordFinder->SetDeltaChord( deltaChord );
//
//    G4double deltaOneStep      = 0.01*mm;
//    fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);
//
//    G4double deltaIntersection = 0.1*mm;
//    fieldManager->SetDeltaIntersection(deltaIntersection);
//
//    G4TransportationManager* transportManager =
//    G4TransportationManager::GetTransportationManager();
//
//    G4PropagatorInField* fieldPropagator =
//    transportManager->GetPropagatorInField();
//
//    G4double epsMin            = 2.5e-5*mm;
//    G4double epsMax            = 0.05*mm;
//
//    fieldPropagator->SetMinimumEpsilonStep(epsMin);
//    fieldPropagator->SetMaximumEpsilonStep(epsMax);
//
//    fieldManager->SetChordFinder(chordFinder);

                                  
//Envelope
   
 G4Box* solidEnv =
   new G4Box("Envelope",                    //its name
        env_sizeXY, env_sizeXY, env_sizeZ); //its size
      
 G4LogicalVolume* logicEnv =
   new G4LogicalVolume(solidEnv,            //its solid
                       world_mat,             //its material
                       "Envelope");
//                       fieldManager);

               
 new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  //  logicEnv->SetFieldManager(fieldManager, true);
  //     
  // Shape 1
  //
  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");
 // G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Si");
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 2*m);
        
  // Cylindrical section shape
  G4double cyl_rad = 1.5*m;
  G4double cyl_height = 2*m;
  G4double shape1_hz = 3.*m;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Tubs* solidShape1 =
    new G4Tubs("Shape1",
     0, cyl_rad, cyl_height, 0, 360*degree);
/*
  // Full sphere shape
  G4double shape1_rmax = 4*cm;
  G4Orb* solidShape1 =    
    new G4Orb("Shape1",                     //its name
              shape1_rmax);                 //its size

  // Sphere shape
  G4double shape1_rmin = 0*cm, shape1_rmax = 4*cm;
  G4double shape1_thetamin = 0.*deg, shape1_thetamax =  180.*deg;    
  G4double shape1_phimin = 0.*deg, shape1_phimax =  360.*deg;    
  G4Sphere* solidShape1 =    
    new G4Sphere("Shape1",                  //its name
        shape1_rmin, shape1_rmax,                //its size
        shape1_phimin, shape1_phimax,            //phi angle
        shape1_thetamin, shape1_thetamax);       //theta angle
     
  // Box shape
  G4double shape1_dx = 8*cm, shape1_dy = 8*cm, shape1_dz = 8*cm;    
  G4Box* solidShape1 =    
    new G4Box("Shape1",                     //its name
         0.5*shape1_dx, 0.5*shape1_dy, 0.5*shape1_dz);     //its size
*/
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
   
    G4double mag = 0.0054;

  //     
  // Shape 2
  //
    a = 60.08*g/mole;
    G4Element*elSO2 = new G4Element(name = "Silicon Dioxide", symbol = "SiO2", z = 12, a);
    
    a = 71.84*g/mole;
    G4Element*elFeOx = new G4Element(name = "Ferrous Oxide", symbol = "FeOx", z = 13, a);
    
    a =40.30*g/mole;
    G4Element*elMgO = new G4Element(name = "Magnesium Oxide", symbol = "MgO", z = 14, a);
    
    a = 56.08*g/mole;
    G4Element*elCaO = new G4Element(name = "Calcium Oxide", symbol = "CaO", z = 15, a);
    
    a = 101.96*g/mole;
    G4Element*elAl2O3 = new G4Element(name = "Aluminum Oxide", symbol = "Al2O3", z = 16, a);
    
    a = 18.02*g/mole;
    G4Element*elH2O = new G4Element(name = "Water", symbol = "H2O", z = 16, a);
    
    density = 1.29*g/cm3;
    
    G4Material* shape2_mat = new G4Material(name = "Soil_Regolith", density, ncomponents = 6);
    shape2_mat->AddElement(elSO2, fractionmass = 52*perCent);
    shape2_mat->AddElement(elFeOx, fractionmass = 21*perCent);
    shape2_mat->AddElement(elMgO, fractionmass = 7*perCent);
    shape2_mat->AddElement(elCaO, fractionmass = 7*perCent);
    shape2_mat->AddElement(elAl2O3, fractionmass = 10*perCent);
    shape2_mat->AddElement(elH2O, fractionmass = 3*perCent);
    
  G4ThreeVector pos2 = G4ThreeVector(0, 0, 10*m);
/*
  //  Shape 2 - conical section shape       
   G4double shape2_rmina =  0.*cm, shape2_rmaxa = 5.*cm;
   G4double shape2_rminb =  0.*cm, shape2_rmaxb = 8.*cm;
   G4double shape2_hz = 3.*cm;
   G4double shape2_phimin = 0.*deg, shape2_phimax = 360.*deg;
   G4Cons* solidShape2 =    
     new G4Cons("Shape2", 
     shape2_rmina, shape2_rmaxa, shape2_rminb, shape2_rmaxb, shape2_hz,
     shape2_phimin, shape2_phimax);
*/

  // Trapezoid shape       
 // G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
//  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
 // G4double shape2_dz  = 6*cm;
 // G4Trd* solidShape2 =
 //   new G4Trd("Shape2",                      //its name
 //             0.5*shape2_dxa, 0.5*shape2_dxb,
 //             0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
    G4double reg2_xy = 16*m;
    G4double reg2_z = 6*m;
    G4Box* solidShape2 =
    new G4Box("Shape2",
              reg2_xy, reg2_xy, reg2_z);
    
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    
    // Hallbach shape
    G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_F");
    G4ThreeVector pos3 = G4ThreeVector(0, 0, -5*m);
    G4double tube_outrad =  1.5*m;
    G4double tube_innerrad = 0.25*m;
    G4double tube_height = 1*m;
    
    //Magnetic Field
    G4UniformMagField* magField =
    new G4UniformMagField(G4ThreeVector(0,mag*tesla,0));
    
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    
    //Set Equation of Motion in Magnetic Field
    
    G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(magField);
    G4MagIntegratorStepper* pStepper = new G4HelixExplicitEuler(fEquation);
    G4MagInt_Driver*pIntgrDriver = new G4MagInt_Driver(1e-5*mm, pStepper, pStepper->GetNumberOfVariables());
    G4ChordFinder* pChordFinder = new G4ChordFinder(pIntgrDriver);
    pChordFinder->SetDeltaChord(1e-3*mm); //Maximum Miss Distance
    fieldMgr->SetChordFinder(pChordFinder);
    fieldMgr->SetDeltaOneStep(1e-3*mm);
    fieldMgr->SetDeltaIntersection(1e-4*mm);
    
    G4Tubs* solidhalbach =
    new G4Tubs("Halbach",                       //its name
               tube_innerrad, tube_outrad, tube_height, 0, 360*degree);//its size
    
    G4LogicalVolume* logichalbach =
    new G4LogicalVolume(solidhalbach,            //its solid
                        shape3_mat,          //its material
                        "Halbach",          //its name
                        fieldMgr);
    
    new G4PVPlacement(0,                       //no rotation
                      pos3,                   //at (0,0,0)
                      logichalbach,                //its logical volume
                      "Halbach",              //its name
                      logicEnv,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    logichalbach->SetFieldManager(fieldMgr, true);
    
  // Set scoring volume to stepping action 
  // (where we will account energy deposit)
  //
  //B1SteppingAction* steppingAction = B1SteppingAction::Instance();
  //steppingAction->SetVolume(logicShape2);
  //steppingAction->SetVolume(logicShape1);
    fScoringVolume = logicShape1;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
