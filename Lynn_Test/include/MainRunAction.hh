//
//  MainRunAction.h
//  
//
//  Created by Lynn Sargent on 1/10/19.
//

#ifndef MainRunAction_h
#define MainRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class MainRunAction : public G4UserRunAction
{
    public:
    MainRunAction();
    virtual ~MainRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
};

#endif
