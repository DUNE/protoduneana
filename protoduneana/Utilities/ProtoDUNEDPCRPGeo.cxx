// ProtoDUNEDPCRPGeo.cxx
//
#include <iostream>
#include <string>

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"

#include "ProtoDUNEDPCRPGeo.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using geo::Point_t;
using protoana::crpgeoinfo;

// ctor
protoana::ProtoDUNEDPCRPGeo::ProtoDUNEDPCRPGeo()
{
  fLemsPerRow = 6;
  fLemWidth   = 50.0; // in cm

  // geometry service
  fGeom = &*art::ServiceHandle<geo::Geometry>();
  fWireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();
}

//
// dtor
protoana::ProtoDUNEDPCRPGeo::~ProtoDUNEDPCRPGeo()
{;}


//
crpgeoinfo protoana::ProtoDUNEDPCRPGeo::GetCRPGeoInfo( Point_t const &pnt ) const
{
  const string myname = "protoana::ProtoDUNEDPCRPGeo::GetCRPGeoInfo: ";

  crpgeoinfo crp_geo_info;

  // find TPC volume
  geo::TPCGeo const* tpc = nullptr;
  try{
    tpc = fGeom->PositionToTPCptr( pnt );
  }
  catch(cet::exception &e){
    cout<<e;
  }

  if( not tpc ) return crp_geo_info;
  crp_geo_info.crpid  = tpc->ID().TPC;

  // determine drift coordinate
  auto const [axis, sign] = tpc->DriftAxisWithSign();
  int tcoord1 = 0;
  int tcoord2 = 0;
  if(axis == geo::Coordinate::X)
    {
      tcoord1 = 1; // Y
      tcoord2 = 2; // Z
    }
  else if(axis == geo::Coordinate::Y)
    {
      tcoord1 = 0; // X
      tcoord2 = 2; // Z
    }
  else // wrong CRP drift
    {
      cout<<myname<<"Bad drift direction "<<axis<<endl;
      return crp_geo_info;
    }


  auto const plane_center = fWireReadoutGeom->Plane({tpc->ID(), 0}).GetCenter();
  crp_geo_info.danode = geo::vect::coord(plane_center, to_int(axis)) - geo::vect::coord(pnt, to_int(axis));

  // point coordinate in the plane
  float planeX = geo::vect::coord(pnt, tcoord2) - geo::vect::coord(plane_center, tcoord2);
  float planeY = geo::vect::coord(pnt, tcoord1) - geo::vect::coord(plane_center, tcoord1);

  // cout<<"Input : "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
  // cout<<"Plane : "
  //     <<tpc->PlaneLocation(0)[0]<<" "
  //     <<tpc->PlaneLocation(0)[1]<<" "
  //     <<tpc->PlaneLocation(0)[2]<<endl;
  // cout<<"Danode : "<<crp_geo_info.danode<<endl;
  // cout<<"Plane coord : "<<planeX<<" "<<planeY<<endl;

  // distance to the CRP edge
  float activeHalfWidth = 0.5 * fLemsPerRow * fLemWidth;
  crp_geo_info.dedge = std::min( (activeHalfWidth - std::abs(planeX)),
                                 (activeHalfWidth - std::abs(planeY)) );
  if( crp_geo_info.dedge < 0 || crp_geo_info.dedge > activeHalfWidth ){
    cout<<myname<<"Bad distance to CRP edge "<<crp_geo_info.dedge<<endl;
    return crp_geo_info;
  }

  // 2d index of each LEM: CRP edges should be about -150.0 -> 150.0 cm
  int LemsHalfRow = fLemsPerRow / 2;
  int icol = int(planeX / fLemWidth); // column index along Z axis
  if( planeX < 0 ) icol += (LemsHalfRow - 1);
  else icol += (LemsHalfRow);

  int irow = int(planeY / fLemWidth); // row index in other coordinate
  if( planeY < 0 ) irow += (LemsHalfRow - 1);
  else irow += (LemsHalfRow);

  crp_geo_info.lemid = icol * fLemsPerRow + irow;
  if( crp_geo_info.lemid >= fLemsPerRow * fLemsPerRow || crp_geo_info.lemid < 0 ){
    // bad LEM index;
    cout<<myname<<"Bad LEM index "<<crp_geo_info.lemid<<endl;
    return crp_geo_info;
  }

  // assumed center of each LEM
  float lemX = (icol + 0.5) * fLemWidth - activeHalfWidth;
  float lemY = (irow + 0.5) * fLemWidth - activeHalfWidth;

  // position wrt to LEM center
  float inlemX = planeX - lemX;
  float inlemY = planeY - lemY;
  crp_geo_info.dlem = std::min( (0.5 * fLemWidth - std::abs(inlemX)),
                                (0.5 * fLemWidth - std::abs(inlemY)) );
  // cout<<"LEM IDs : "<<icol<<" "<<irow<<" "<<crp_geo_info.lemid<<"\n";
  // cout<<"LEM pos : "<<lemX<<" "<<lemY<<"\n";
  // cout<<"In LEM pos :  "<<inlemX<<" "<<inlemY<<"\n";
  // cout<<"Distance to LEM border : "<<crp_geo_info.dlem<<endl;

  if( crp_geo_info.dlem < 0 || crp_geo_info.dlem > 0.5*fLemWidth ){
    cout<<myname<<"Bad distance to LEM border "<<crp_geo_info.dlem<<endl;
    return crp_geo_info;
  }

  crp_geo_info.valid = true;

  return crp_geo_info;
}
