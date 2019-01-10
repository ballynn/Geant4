//
//  MainDetectorConstruction.h
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#ifndef MainDetectorConstruction_h
#define MainDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4Element;
class G4Isotope;

/// Detector construction class to define materials and geometry.

class MainDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
    MainDetectorConstruction();
    virtual ~MainDetectorConstruction();
    
    public:
    virtual G4VPhysicalVolume* Construct();
};

#endif /* MainDetectorConstruction_h */
