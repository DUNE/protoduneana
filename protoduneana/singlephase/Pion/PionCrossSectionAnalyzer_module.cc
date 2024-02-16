////////////////////////////////////////////////////////////////////////
// Class:       PionCrossSectionAnalyzer
// Plugin Type: analyzer (art v3_01_02)
// File:        PionCrossSectionAnalyzer_module.cc
//
// Generated at Wed Mar  6 14:10:55 2019 by Jacob Calcutt using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
//#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
//#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
//#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
//#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

#include "art_root_io/TFileService.h"

// ROOT includes
#include "TTree.h"
#include "TH1D.h"

class PionCrossSectionAnalyzer;


class PionCrossSectionAnalyzer : public art::EDAnalyzer {
public:
  explicit PionCrossSectionAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionCrossSectionAnalyzer(PionCrossSectionAnalyzer const&) = delete;
  PionCrossSectionAnalyzer(PionCrossSectionAnalyzer&&) = delete;
  PionCrossSectionAnalyzer& operator=(PionCrossSectionAnalyzer const&) = delete;
  PionCrossSectionAnalyzer& operator=(PionCrossSectionAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;
  
  double distance(double x1, double y1, double z1, double x2, double y2, double z2);

  void reset();

private:

  TH1D *hInteractingKEInel; 
  TH1D *hInteractingKEAbs;
  TH1D *hInteractingKECex;
  TH1D *hInteractingKEQE;
  TH1D *hInteractingKEDCex;
  TH1D *hInteractingKEProd;
  TH1D *hIncidentKE;

  double minX =  -360.0;
  double maxX = 360.0;
  double minY =0.0;
  double maxY = 600.0;
  double minZ =  0.0; // G10 at z=1.2
  double maxZ = 695.0;

  std::string fGeneratorTag;

};

double PionCrossSectionAnalyzer::distance(double x1, double y1, double z1, double x2, double y2, double z2){
  double d = TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return d;
}

PionCrossSectionAnalyzer::PionCrossSectionAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}   ,
  fGeneratorTag(p.get<std::string>("GeneratorTag"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void PionCrossSectionAnalyzer::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  hIncidentKE        = tfs->make<TH1D>("hIncidentKE"   , "Incident Kinetic Energy [MeV]"   , 420, -100, 2000); 
  hInteractingKEInel = tfs->make<TH1D>("hInteractingKEInel", "Inelastic Interacting Kinetic Energy [MeV]", 420, -100, 2000);
  hInteractingKEQE = tfs->make<TH1D>("hInteractingKEQE", "Quasi-elastic Interacting Kinetic Energy [MeV]", 420, -100, 2000);
  hInteractingKEAbs = tfs->make<TH1D>("hInteractingKEAbs", "Absorption Interacting Kinetic Energy [MeV]", 420, -100, 2000);
  hInteractingKECex = tfs->make<TH1D>("hInteractingKECex", "Charge Exchange  Interacting Kinetic Energy [MeV]", 420, -100, 2000);
  hInteractingKEDCex = tfs->make<TH1D>("hInteractingKEDCex", "Double Charge Exchange  Interacting Kinetic Energy [MeV]", 420, -100, 2000);
  hInteractingKEProd = tfs->make<TH1D>("hInteractingKEProd", "Pion Production  Interacting Kinetic Energy [MeV]", 420, -100, 2000);

}
void PionCrossSectionAnalyzer::endJob(){

}

void PionCrossSectionAnalyzer::analyze(art::Event const& evt){
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList(); 

  // Get the truth utility to help us out
  protoana::ProtoDUNETruthUtils truthUtil;
  // Firstly we need to get the list of MCTruth objects from the generator. The standard protoDUNE
  // simulation has fGeneratorTag = "generator"
  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
  // mcTruths is basically a pointer to an std::vector of simb::MCTruth objects. There should only be one
  // of these, so we pass the first element into the function to get the good particle
  const simb::MCParticle* geantGoodParticle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

  bool keepInteraction = false;

  if(geantGoodParticle != 0x0){
    std::cout << "Found GEANT particle corresponding to the good particle with pdg = " << geantGoodParticle->PdgCode() << std::endl;
    std::cout << "ID: " << geantGoodParticle->TrackId() << std::endl;
    std::cout << "Mother: " << geantGoodParticle->Mother() << std::endl;

    if(geantGoodParticle->PdgCode()==211 && geantGoodParticle->Process()=="primary" && geantGoodParticle->NumberDaughters()!=0){
      std::cout<<"Track ID of the geantGoodParticle "<<geantGoodParticle->TrackId()<<std::endl;
      std::cout<<"Process:"<<geantGoodParticle->Process()<<std::endl;
      
      simb::MCTrajectory truetraj=geantGoodParticle->Trajectory();

      auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
      geo::View_t view = wireReadout.Plane({0, 0, 2}).View();
      auto simIDE_prim=bt_serv->TrackIdToSimIDEs_Ps(geantGoodParticle->TrackId(),view);
      std::map<double, sim::IDE> orderedSimIDE;
      for (auto& ide : simIDE_prim) orderedSimIDE[ide->z]= *ide;

      std::string interactionLabel="";
      double mass=geantGoodParticle->Mass();
      //Store the kinetic energy and momentum on z at WC4. Just for cross check 
      auto inTPCPoint  = truetraj.begin(); 

      //--------------------------------------------------------
      // Identify the first trajectory point inside the TPC
      // Loop From First TrajPoint --> First Point in TPC 
      // Stop when you get into the TPC
      for ( auto t = truetraj.begin(); t!= std::prev(truetraj.end()); t++){
        auto pos = t->first;
        if (pos.Z() < minZ) continue;
        else if (pos.X() < minX || pos.X() > maxX ) continue;
        else if (pos.Y() < minY || pos.Y() > maxY ) continue;
        else {
          inTPCPoint = t;
          break;
        }
      }

      if (inTPCPoint !=truetraj.begin()) {

        // Identify the last interesting trajectory point in TPC
        auto finTPCPoint = std::prev(truetraj.end()); 
        // The last point is a bit more complicated:
        // if there's no interaction, then it is simply the last point in the TPC
        // if there's one or more interaction points, it's the first interaction point deemed interesting (no coulomb)
        // Take the interaction Map... check if there's something there
        auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
        std::cout<<"thisTrajectoryProcessMap size: "<<thisTracjectoryProcessMap.size()<<std::endl;
    
        int inel=0;
        for(auto const& couple: thisTracjectoryProcessMap) {
          // I'm not interested in the CoulombScattering, discard this case
          if ((truetraj.KeyToProcess(couple.second)).find("CoulombScat")!= std::string::npos) continue;
          
          // Let's check if the interaction is in the the TPC
          auto     interactionPos4D =  (truetraj.at(couple.first)).first ;        
          if      (interactionPos4D.Z() <  minZ || interactionPos4D.Z() > maxZ ) continue;
          else if (interactionPos4D.X() <  minX || interactionPos4D.X() > maxX ) continue;
          else if (interactionPos4D.Y() <  minY || interactionPos4D.Y() > maxY ) continue;
          if ((truetraj.KeyToProcess(couple.second)).find("InElastic")!= std::string::npos) inel++;

          // If we made it here, then this is the first interesting interaction in the TPC
          // Our job is done!!! Great! Store the interaction label and the iterator for the final point
          std::cout << "key: " << size_t(couple.second) << std::endl;
          interactionLabel = truetraj.KeyToProcess(couple.second);
          std::cout<<"interaction Label "<<interactionLabel<<" EndProcess "<<geantGoodParticle->EndProcess()<<std::endl;
          finTPCPoint = truetraj.begin() + couple.first; 
          keepInteraction=true;
          break; 
        }

        sim::ParticleList::const_iterator itPart=plist.begin();
  
        // If I didn't find anything interesting in the intereaction map, let's loop back!
        if (!keepInteraction) {
          //Loop on the daughters 
          for(size_t iPart1=0;(iPart1<plist.size()) && (plist.begin()!=plist.end());++iPart1){
            const simb::MCParticle* pPart=(itPart++)->second;
            //   art::Ptr<simb::MCTruth> const& mcDaught=pi_serv->ParticleToMCTruth_P(pPart)
            //We keep only the dauthers of the primary not coming from elastic or inelastic scattering
            if (pPart->Mother()  != 1 ) continue;
            //  std::cout<<"pPart mother "<<pPart->Mother()<<std::endl;
            if ((pPart->Process()).find("astic")!= std::string::npos) continue;
            if ((pPart->Process()).find("CoulombScat")!= std::string::npos) continue;
            //Is the daughter born inside the TPC? If yes, store the process which created it 
            simb::MCTrajectory trueDaugthTraj = pPart->Trajectory();              
            if (trueDaugthTraj.begin()->first.Z() < minZ || trueDaugthTraj.begin()->first.Z() > maxZ) continue;
            else if (trueDaugthTraj.begin()->first.X() <   minX || trueDaugthTraj.begin()->first.X() > maxX ) continue;
            else if (trueDaugthTraj.begin()->first.Y() <   minY || trueDaugthTraj.begin()->first.Y() > maxY ) continue;
            else {
              interactionLabel = pPart->Process();
          
              break;
            }
          }
          for ( auto t = std::prev(truetraj.end()); t!= truetraj.begin(); t--){
            auto pos = t->first;
          
            if (pos.Z() > maxZ) continue;
            else if (pos.X() <   minX || pos.X() > maxX ) continue;
            else if (pos.Y() <   minY || pos.Y() > maxY ) continue;
            else {
              finTPCPoint = t;
              break;
            }
          }
        }

        if (finTPCPoint != inTPCPoint){
          auto posFin = finTPCPoint->first;
          auto posIni = inTPCPoint->first;
          //Let's record what the initial and final points are.

          auto totLength = distance(posFin.X(), posFin.Y(), posFin.Z(),posIni.X(), posIni.Y(), posIni.Z() );

    
          // We want to chop up the points between the first and last uniformly
          // and ordered by Z
          // Order them in order of increasing Z
          std::map<double, TVector3> orderedUniformTrjPts;
          // We want the first and uniform point to coincide with the 
          // the first and last points we just found 
          auto positionVector0 = (inTPCPoint ->first).Vect(); 
          auto positionVector1 = (finTPCPoint->first).Vect(); 
          orderedUniformTrjPts[positionVector0.Z()] = positionVector0;
          orderedUniformTrjPts[positionVector1.Z()] = positionVector1;

          const double trackPitch = 0.4792;
          // I do have space for at least one extra point, so let's put it there!
          // Calculate how many extra points I need to put between the new first point and the second TrajPoint
          int    nPts            = (int) (totLength/trackPitch);
          for (int iPt = 1; iPt <= nPts; iPt++ ){
            auto newPoint = positionVector0 + iPt*(trackPitch/totLength) * (positionVector1 - positionVector0);
            orderedUniformTrjPts[newPoint.Z()] = newPoint;
          }

          // If the distance between the last point and the second to last is less then 0.235
          // eliminate the second to last point
          auto lastPt         = (orderedUniformTrjPts.rbegin())->second;
          auto secondtoLastPt = (std::next(orderedUniformTrjPts.rbegin()))->second;
          double lastDist = distance(lastPt.X(),lastPt.Y(),lastPt.Z(),secondtoLastPt.X(),secondtoLastPt.Y(),secondtoLastPt.Z());
          if (lastDist < 0.240){
            orderedUniformTrjPts.erase((std::next(orderedUniformTrjPts.rbegin()))->first );
          } 
    
          // Calculate the initial kinetic energy
          auto initialMom =     inTPCPoint->second;
          double initialKE = 1000*(TMath::Sqrt(initialMom.X()*initialMom.X() + initialMom.Y()*initialMom.Y() + initialMom.Z()*initialMom.Z() + mass*mass ) - mass); 
          double kineticEnergy = initialKE;
          auto old_it = orderedUniformTrjPts.begin();
          for (auto it = std::next(orderedUniformTrjPts.begin()); it != orderedUniformTrjPts.end(); it++, old_it++ ){
          
            auto oldPos        = old_it->second;
            auto currentPos    =     it->second;
        
            double uniformDist =  (currentPos - oldPos).Mag();
          
            //Calculate the energy deposited in this slice          
            auto old_iter = orderedSimIDE.begin();
            double currentDepEnergy = 0.;
            for ( auto iter= orderedSimIDE.begin(); iter!= orderedSimIDE.end(); iter++,old_iter++){
              auto currentIde = iter->second;
              // std::cout<<"Z position of the trajectory hits "<<currentIde.z<<std::endl;
              if ( currentIde.z < oldPos.Z()) continue;
              if ( currentIde.z > currentPos.Z()) continue;
              currentDepEnergy += currentIde.energy;
            }

            // avoid overfilling super tiny energy depositions
            if (currentDepEnergy/uniformDist < 0.1 ) continue;

            //Calculate the current kinetic energy
            kineticEnergy -= currentDepEnergy;
            hIncidentKE->Fill(kineticEnergy);
          }

          if(interactionLabel.find("Inelastic")!= std::string::npos ){
            hInteractingKEInel->Fill(kineticEnergy);

            int nPiPlus_truth = 0;
            int nPiMinus_truth = 0;
            int nPi0_truth = 0;
            int nProton_truth = 0;
            int nNeutron_truth = 0;


            std::cout << "NDaughters: " << geantGoodParticle->NumberDaughters() << std::endl;
            for( int i = 0; i < geantGoodParticle->NumberDaughters(); ++i ){
              std::cout << "Checking " << i << std::endl;
              int daughterID = geantGoodParticle->Daughter(i);

              std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
              //Skip photons, neutrons, the nucleus
              if( plist[ daughterID ]->PdgCode() == 22 || /*plist[ daughterID ]->PdgCode() == 2112 ||*/ plist[ daughterID ]->PdgCode() > 1000000000  ) continue;

              auto part = plist[ daughterID ];
              int pid = part->PdgCode();
              std::cout << "PID: " << pid << std::endl;
              std::cout << "Process: " << part->Process() << std::endl;

              if( pid == 211 )  nPiPlus_truth++;
              else if( pid == -211 ) nPiMinus_truth++;
              else if( pid == 111 )  nPi0_truth++;
              else if( pid == 2212 ) nProton_truth++;
              else if( pid == 2112 ) nNeutron_truth++;

            }

            if( nPiPlus_truth == 1 && nPiMinus_truth == 0 && nPi0_truth == 0 ){
              std::cout << "Filling QE" << std::endl;
              hInteractingKEQE->Fill(kineticEnergy); 
            }
            else if( nPiPlus_truth == 0 && nPiMinus_truth == 0 && nPi0_truth == 0 ){
              std::cout << "Filling Abs" << std::endl;
              hInteractingKEAbs->Fill(kineticEnergy); 
            }
            else if( nPiPlus_truth == 0 && nPiMinus_truth == 0 && nPi0_truth == 1 ){
              std::cout << "Filling Cex" << std::endl;
              hInteractingKECex->Fill(kineticEnergy); 
            }
            else if( nPiPlus_truth == 0 && nPiMinus_truth == 1 && nPi0_truth == 0 ){
              std::cout << "Filling DCex" << std::endl;
              hInteractingKEDCex->Fill(kineticEnergy); 
            }
            else{
              std::cout << "Filling Prod" << std::endl;
              hInteractingKEProd->Fill(kineticEnergy); 
            }
          }
    
          keepInteraction = false;

          }
        }
      }



  }

}

DEFINE_ART_MODULE(PionCrossSectionAnalyzer)
