//
//  MainDetectorConstruction.c
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#include "MainDetectorConstruction.hh"
#include "MainSteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

// I believe this initializes the detector construction with Geant4 Library
MainDetectorConstruction::MainDetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//I'm not sure what this is doing
MainDetectorConstruction::~MainDetectorConstruction()
{ }

//Construct World
G4VPhysicalVolume* MainDetectorConstruction::Construct()
{

    //
    //Create World
    //
    G4double world_sizeXY = 10.0*m;
    G4double world_sizeZ = 10.0*m;

    //Define Mars Atmosphere
    G4double z, a, fractionmass, density, pressure;
    G4String name, symbol;
    G4int ncomponents;
    
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
    
    density = 3.93*g/cm3;
    pressure = 840*kg/(m*s);                                                        //reference A. Firan Paper
    G4Material* world_mat = new G4Material(name = "Mars_Air", density, pressure, ncomponents = 5);
    world_mat->AddElement(elN, fractionmass = 1.89*perCent);
    world_mat->AddElement(elCO2, fractionmass = 96.0*perCent);
    world_mat->AddElement(elAr, fractionmass = 1.93*perCent);
    world_mat->AddElement(elO, fractionmass = 0.145*perCent);
    world_mat->AddElement(elCO, fractionmass = 0.1*perCent);

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;

    G4Box* solidWorld =
        new G4Box("World",                                                  //name
                  0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //size of world

    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld,                                     //its solid
                            world_mat,                                      //its material
                            "World");                                       //its name
    G4VPhysicalVolume* physWorld =
        new G4PVPlacement(0,                                                //no rotation
                          G4ThreeVector(),                                  //at (0,0,0)
                          logicWorld,                                       //its logical volume
                          "World",                                          //its name
                          0,                                                //its mother volume
                          false,                                            //no boolean operation
                          0,                                                //copy number
                          checkOverlaps);                                   //overlaps checking

    //
    //CylindricalEnvelope
    //
    G4double cyl_rad = 1*m;
    G4double cyl_height = 2*m;

    G4Tubs*env_cyl =
        new G4Tubs("Envelope_Cylinder",                                     //name
                   0, cyl_rad, cyl_height, 0, 360*degree);                  //size

    G4LogicalVolume* logicEnv =
        new G4LogicalVolume(env_cyl,                                        //its solid
                            world_mat,                                      //its material
                            "Envelope_Cylinder");                           //its name
    
        new G4PVPlacement(0,                                                //no rotation
                          G4ThreeVector(),                                  //at 0,0,0
                          logicEnv,                                         //its logical volume
                          "Envelope_Cylinder",                              //its name
                          logicWorld,                                       //its mother volume
                          false,                                            //no boolean operation
                          0,                                                //copy number
                          checkOverlaps);                                   //overlap checking
    
    //
    // Mars Regolith - loose soil model (See Ana Firan Presentation Boulder)
    //
    G4double reg1_xy = 0.5*world_sizeXY;
    G4double reg1_z = 0.5*cm;
    G4ThreeVector pos1 = G4ThreeVector(0, 0, -0.5*cyl_height);
    
    G4Box* Soil_Mars_Reg1 =
        new G4Box("Mars_Regolith1",
                  reg1_xy, reg1_xy, reg1_z);

    //Define Mars Regolith Material
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
    
    G4Material*Mars_Reg1 =
        new G4Material(name = "Soil_Regolith", density, ncomponents = 6);
    Mars_Reg1->AddElement(elSO2, fractionmass = 52*perCent);
    Mars_Reg1->AddElement(elFeOx, fractionmass = 21*perCent);
    Mars_Reg1->AddElement(elMgO, fractionmass = 7*perCent);
    Mars_Reg1->AddElement(elCaO, fractionmass = 7*perCent);
    Mars_Reg1->AddElement(elAl2O3, fractionmass = 10*perCent);
    Mars_Reg1->AddElement(elH2O, fractionmass = 3*perCent);

    G4LogicalVolume*Soil_Mars_Reg =
        new G4LogicalVolume(Soil_Mars_Reg1,
                            Mars_Reg1,
                            "Soil_Mars_Regolith");
    
        new G4PVPlacement(0,                                                //no rotation
                          pos1,                                             //at 0,0,0
                          logicEnv,                                         //its logical volume
                          "Soil_Mars_Regolith",                              //its name
                          logicWorld,                                       //its mother volume
                          false,                                            //no boolean operation
                          0,                                                //copy number
                          checkOverlaps);                                   //overlap checking)
    
    //
    // Mars Regolith - rock soil model (See Ana Firan Presentation Boulder)
    //
    G4double reg2_xy = 0.5*world_sizeXY;
    G4double reg2_z = 0.5*world_sizeZ - 0.5*cyl_height - reg1_z;
    G4ThreeVector pos2 = G4ThreeVector(0, 0, -0.5*cyl_height - reg1_z);
    
    G4Box* Solid_Mars_Reg2 =
    new G4Box("Mars_Regolith2",
              reg2_xy, reg2_xy, reg2_z);
    
    //Rock Material
    density = 3.49*g/cm3;
    G4Material*Mars_Reg2 =
        new G4Material(name = "Rock_Regolith", density, ncomponents = 6);
    Mars_Reg2->AddElement(elSO2, fractionmass = 52*perCent);
    Mars_Reg2->AddElement(elFeOx, fractionmass = 21*perCent);
    Mars_Reg2->AddElement(elMgO, fractionmass = 7*perCent);
    Mars_Reg2->AddElement(elCaO, fractionmass = 7*perCent);
    Mars_Reg2->AddElement(elAl2O3, fractionmass = 10*perCent);
    Mars_Reg2->AddElement(elH2O, fractionmass = 3*perCent);

    G4LogicalVolume*Solid_Mars_Reg =
        new G4LogicalVolume(Solid_Mars_Reg2,
                            Mars_Reg2,
                            "Rock_Mars_Regolith");
    
        new G4PVPlacement(0,                                                //no rotation
                          pos2,                                             //at pos2
                          logicEnv,                                         //its logical volume
                          "Rock_Mars_Regolith",                              //its name
                          logicWorld,                                       //its mother volume
                          false,                                            //no boolean operation
                          0,                                                //copy number
                          checkOverlaps);                                   //overlap checking)
    
    // Set scoring volume to stepping action
    // (where we will account energy deposit)
    //
    MainSteppingAction* steppingAction = MainSteppingAction::Instance();
    ////steppingAction->SetVolume(logicShape1);
    steppingAction->SetVolume(logicWorld);
    
    
    //
    //always return the physical World
    //
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
