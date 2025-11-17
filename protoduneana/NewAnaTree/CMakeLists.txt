#art_make(BASENAME_ONLY
#  MODULE_LIBRARIES
#  cetlib::cetlib
#  cetlib_except::cetlib_except
#  larsim::Simulation
#  lardataobj::Simulation
#  lardataobj::RecoBase
#  lardata::Utilities
#  lardata::ArtDataHelper
#  larcorealg::Geometry
#  nusimdata::SimulationBase
#  larsim::MCCheater_PhotonBackTrackerService_service           
#  larsim::MCCheater_BackTrackerService_service           
#  larsim::MCCheater_ParticleInventoryService_service
#  larreco::RecoAlg
#  ROOT::Core ROOT::Tree ROOT::Hist
#  ProtoDUNEUtilities
#  geant4reweight::ReweightBaseLib
#  geant4reweight::PropBaseLib
#  art::Framework_Core
#  art::Framework_Principal
#  art::Framework_Services_Registry
#  art_root_io::tfile_support
#  ROOT::Core
#  art_root_io::TFileService_service
#  art::Persistency_Common
#  canvas::canvas
#  art::Persistency_Provenance
#  art::Utilities
#  messagefacility::MF_MessageLogger
#  hep_hpc_hdf5
#  hep_hpc_Utilities
#  ${HDF5_LIBRARIES}
#  )
include_directories($ENV{SRPROXY_INC})
include_directories($ENV{DUNEANAOBJ_INC})
cet_build_plugin(CAFMakerPDUNE   art::module LIBRARIES
                        duneanaobj::StandardRecord
                        duneanaobj::StandardRecordFlat
                        nugen::NuReweight
                        nugen::NuReweight_art
			nugen::EventGeneratorBase_GENIE
			nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Core ROOT::Hist ROOT::Tree
			${GENIE_LIB_LIST}
                        systematicstools::interface
			systematicstools::interpreters
			systematicstools::utility
                        dunereco::CVN_func
			larpandora::LArPandoraInterface
                        larreco::Calorimetry
                        ProtoDUNEUtilities
                        dunereco::AnaUtils
                        larsim::Utils
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
