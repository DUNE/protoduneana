// #include "gallery/Event.h"
// #include "canvas/Utilities/InputTag.h"
// #include "nusimdata/SimulationBase/MCParticle.h"

// void readmc() {
//   std::vector<std::string> files{"pdhd_prod_beam__gen_g4_IonScintPDExt_PDInt_266463_7_1_20251211T023321Z_detsim_reco1.root"};
//   art::InputTag tag("largeant", "", "G4Stage1");
 
//   TH1F *hm = new TH1F("hm", ";Mother track ID;Particles", 1001, -0.5, 1000.5);
//   // TH1F *h1 = new TH1F("h1", ";Track ID;Particles", 1001, -0.5, 1000.5); 
//   TH1F *hid = new TH1F("hid", ";Track ID (mother = 0);Particles", 10001, -0.5, 10000.5);
//   TH1F *hz = new TH1F("hz", ";Start z [cm];Particles", 150, -3000, 3000);
//   TH1F *hy = new TH1F("hy", ";Start y [cm];Particles", 150, -3000, 3000);
//   TH1F *hx = new TH1F("hx", ";Start x [cm];Particles", 150, -3000, 3000);
//   TH1F *hpdg = new TH1F("hpdg", ";PDG code;Particles", 100, -500, 2500);
//   //TH1F *hzBkg  = new TH1F("hzBkg",  ";Start z [cm];Particles", 150, -2000, 2000);
  
//   const int targetMother = 0;
//   const int targetID = 100;
//   std::vector<int> nPerEvent;
//   std::vector<int> mPerEvent;
//   std::vector<int> TrackId;
//   int iEvt = 0;
 
//   for (gallery::Event ev(files); !ev.atEnd(); ev.next(), ++iEvt) {
//     auto const& parts = *ev.getValidHandle<std::vector<simb::MCParticle>>(tag);
//     if (iEvt == 0) {
//       int shown = 0;
//       for (auto const& p : parts) {
//         std::cout << p << std::endl;
//         if (++shown == 10) break;
//       }
//     }

//     int n = 0;
//     int m = 0;

//     for (auto const& p : parts) {
//       hm->Fill(p.Mother());
//       // h1->Fill(p.TrackId());
//       if (p.Mother() == targetMother) {
//         hid->Fill(p.TrackId());
//         hz->Fill(p.Vz(0));
//         hy->Fill(p.Vy(0));
//         hx->Fill(p.Vx(0));
//         hpdg->Fill(p.PdgCode());
//         TrackId.push_back(p.TrackId());
//         ++n;
//         if (p.TrackId() == targetID) {
// 	        ++m;
//         }
//         // if (p.TrackId() == 1 && std::abs(p.PdgCode()) == 211) hzBeam->Fill(p.Vz(0));
//         // else if (p.TrackId() == 1) hzBkg->Fill(p.Vz(0));
//         }
//       }
//     nPerEvent.push_back(n);
//     mPerEvent.push_back(m);
//     std::cout << "Event " << iEvt << ": " << n
//               << " particles with mother " << targetMother << std::endl;
//     std::cout << "Event " << iEvt << ": " << m
//               << " particles with mother " << targetMother << " and " << targetID << std::endl;
//   }
//   std::cout << "Total objects stored: " << TrackId.size() << std::endl;
//   // for (int v : mPerEvent) std::cout << v << " ";
//   // std::cout << std::endl;

//   hm->GetXaxis()->SetRangeUser(0, 10);
//   //hpdg->GetXaxis()->SetRangeUser(-40, 40);
//   new TCanvas("cm"); hm->Draw("HIST");
//   // new TCanvas("c2"); h1->Draw("HIST");
//   new TCanvas("cid"); hid->Draw("HIST");
//   new TCanvas("cz"); hz->Draw("HIST");
//   new TCanvas("cy"); hy->Draw("HIST");
//   new TCanvas("cx"); hx->Draw("HIST");
//   new TCanvas("cpdg"); hpdg->Draw("HIST");
//   //new TCanvas("czbkg");  hz->Draw("HIST");
// }
