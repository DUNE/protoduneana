#ifndef PROTODUNEDP_CRPGEO_H
#define PROTODUNEDP_CRPGEO_H

///////////////////////////////////////////////////////////////
//
// Utility functions to access geo uniformation for 
// ProtoDUNE CRPs
//
// The positions of each CRP are obtained from the
// geometry service
//
// For the LEM numbering we follow the convention
//
//  |5 |11|   |35|
//  |4 |10|   |34|
//  |3 |9 |   |33|
//  |2 |8 |   |32|
//  |1 |7 |   |31|
//  |0 |6 |...|30|
//  |--------------------> Z axis
//  
// Created: vglaymov@ipnl.in2p3.fr
//
///////////////////////////////////////////////////////////////

//#include "fhiclcpp/ParameterSet.h"
#include "larcorealg/Geometry/fwd.h"

#include <utility>
//#include <vector>

namespace protoana {
  
  typedef struct crpgeoinfo {
    bool  valid;  // flag
    int   crpid;  // CRP ID == TPC ID
    int   lemid;  // LEM ID
    float danode; // distance to CRP along drift
    float dlem;   // min distance to LEM border
    float dedge;  // min distance to CRP edge border
  crpgeoinfo() : valid(false), crpid(-1), lemid(-1), danode(0), dlem(0), dedge(0)
    {;}
  } crpgeoinfo;
  
  class ProtoDUNEDPCRPGeo
    {
    
    public:
      ProtoDUNEDPCRPGeo();
      ~ProtoDUNEDPCRPGeo();
    
      crpgeoinfo GetCRPGeoInfo( geo::Point_t const &pnt) const;
      int GetLemsPerCRP() const 
      { return fLemsPerRow * fLemsPerRow; } 

    private:
      //int fLogLevel;
      int fLemsPerRow; //
      float fLemWidth; // cm
      geo::GeometryCore const* fGeom;
      geo::WireReadoutGeom const* fWireReadoutGeom;
    };



} // namespace

#endif
  

