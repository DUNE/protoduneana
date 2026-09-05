#Created by Shu, 20260129---

Example of command:
lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251101/cosmic_reco_run101_job246_eventID1231_20251225T034125Z_19482.root -e 101:0:1231 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251102/cosmic_reco_run102_job1976_eventID9881_20251230T065540Z_15759.root -e 102:0:9882 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251101/cosmic_reco_run101_job485_eventID2426_20251225T043137Z_20013.root -e 101:0:2427 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251102/cosmic_reco_run102_job696_eventID3481_20251230T041123Z_4272.root -e 102:0:3481 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251101/cosmic_reco_run101_job913_eventID4566_20251225T052801Z_9884.root -e 101:0:4568 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251101/cosmic_reco_run101_job85_eventID426_20251225T031346Z_27758.root -e 101:0:429 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251101/cosmic_reco_run101_job246_eventID1231_20251225T034125Z_19482.root -e 101:0:1232 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251102/cosmic_reco_run102_job1592_eventID7961_20251230T060613Z_22808.root -e 102:0:7965 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251102/cosmic_reco_run102_job522_eventID2611_20251230T035549Z_3227.root -e 102:0:2614 -n 1

lar -c peakBackTrack.fcl -s /pnfs/dune/scratch/users/szh2/MC_pdhd_Michel/artFile_20251102/cosmic_reco_run102_job1802_eventID9011_20251230T063349Z_10622.root -e 102:0:9013 -n 1


lar -c peakBackTrack.fcl -s /exp/dune/data/users/szh2/running_results/MC_PDHD_list/artFiles/official_Production2025/artFile_examples/pdhd_prod_beam__gen_g4_IonScintPDExt_PDInt_259388_167_1_20251209T091115Z_detsim_reco1.root -e 20250627:1:479 -n 1



#Commands to get runID:subRunID:eventID:
#Shu, 20260818--

1.  root -l <root file>
2.  [inside ROOT]: .ls
3.  [inside ROOT]: Events->Scan("EventAuxiliary.id_.subRun_.run_.run_:EventAuxiliary.id_.subRun_.subRun_:EventAuxiliary.id_.event_")

