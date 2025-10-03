# SingleHit Analysis Module

This module (`SingleHit_module.cc`) is part of theDUNE/protoduneana/verticaldrift/singlehit fork (same code than duneana/duneana/CalibAna).  
It processes reconstructed TPC hits from given generator (gaushit ….), reconstruct hits position (creation of spacepoints) via double or triple time coincidence between triggered wires, identifies isolated (Z,Y, time) hits, clusters possible decays, and stores results in two ROOT trees.

---

##  Main Functions

- **`beginJob()`** 
	– Initialize ROOT trees and branches.  
	- Initialise working variable.
	- Retreive information form configuration files (fcl file)
	- Retreive geometric informations and detector informations (clock, dimension)

- **`analyze(const art::Event &event)`** – Main event loop:  
	- Reads hits (and loop on all collection hits) from the event:
		- retreive charge, time etc 
		- if MC retreive the simulation information (generator, true energy …)
		- With time information of a given hit (triggered collection wire) we form the list of wires (on the two other plane) that are in time coincidence (triggered in a given window given in the fcl file)
		- Then using geometrical function we check that the wire in time coincidence really cross the original collection wire from which the hit originate, from that we retreive (for every collection hit) a list of crossing wire on the two other planes
		- At this point, we can create a list of Spacepoint for each hits (position of crossing point)
        - The first output TTree (hitTree) is filled at this moment with 1 entry = 1 hit, this tree is not for the analysis but more for effecting checks

	- We start the cluster analysis, we fill lot of working vector (*_byEvents). These vectors store the informations of all spacepoints within an event, with all information being at the i-th component for the i-th spacepoint
	-  Now we create the vector vIso, the vector of indices of isolated spacepoints (applying the double sphere veto)
	- From this vIso vector, we start clustering points using app algorithm
	- Then the hit information are retreive within one cluster to define a position (energy weighted barycentre) and energy (sum of hits) on each plane
	- The Clustertree is then filled, one entry = vector of cluster for event 

- **`endJob()`** – Finalization (if needed).  

---

##  Output Trees

The module produces **two ROOT trees**:

1. **`HitTree`** – Contains isolated hit information:  
   - Run, Subrun, Event  
   - Plane, Wire  
   - PeakTime, RMS  
   - Charge, Amplitude  
   - Isolation flag  
   - associated spacepoints

2. **`ClusterTree`** – Contains cluster-level variables:  
   - Cluster size (# hits)  
   - Cluster charge  
   - Cluster position 
   - Associated MC truth (if available)  
   - Generator tags  

##  Helper Functions

The module defines a large number of helper routines. They can be grouped as follows:

###  MC Truth & Utilities
- `GetFirstMCTruthEnergy(...)`  
- `GetGeneratorTag(...)`  

###  Hit Processing @ spacepoint building
- `GetSingle(...)`  
- `NearOrFar(...)`  
- `GetTimeIsolation(...)`  
- `GetListOfTimeCoincidenceHit(...)`  
-  GetListOfCrossingChannel(…)`

###  Crossing & 3D Point Building & isolation
- `GetListOfCrossingChannel(...)`  
- `GetListOf3ViewsPoint(...)`  
- `GetXYZIsolatedPoint(...)`  

###  Clustering & Kpp
- `gen_yz(...)`  
- `gen_yzt(...)`  
- `dist2(...)`  
- `randf(...)`  
- `nearest(...)`  
- `reallocate(...)`  
- `mean(...)`  
- `kpp(...)`  
- `lloyd(...)`  
- `GetData(...)`  
- `CheckCompletude(...)`  
- `CheckClusters(...)`  
- `GetCluster(...)`  

---

### Running

Run the module with:

```bash
Lar -c run_*.fcl -n * reco_file.root -t output_singlehit.root

###
