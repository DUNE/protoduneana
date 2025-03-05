art_make( MODULE_LIBRARIES
                        lardataobj::RawData
			larcorealg::Geometry
                        larcorealg::CoreUtils
			larcore::Geometry_Geometry_service
                        lardataobj::Simulation
                        nusimdata::SimulationBase
			lardataobj::RecoBase
			lardataobj::AnalysisBase
			larsim::MCCheater_ParticleInventoryService_service
			larsim::MCCheater_BackTrackerService_service
			larsim::MCCheater_PhotonBackTrackerService_service
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib 
                        cetlib_except::cetlib_except
			ROOT::Core ROOT::Tree ROOT::Hist
                        ROOT::GenVector
                        ROOT::Geom
                        ProtoDUNEUtilities
)

install_headers()
install_fhicl()
install_source()
install_scripts()
