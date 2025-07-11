#include "qvector_writter.h"
#include "TFile.h"
#include "TTree.h"
#include "HepMC3/Data/GenEventData.h"
#include "smash/pdgcode.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " input.root output.root config.yaml\n";
        return 1;
    }
    double total = 0.;
    double dn_deta = 0.;
    double mean_pt = 0.;
    Config cfg = read_config(argv[3]);
    qvector_writter analyzer(cfg);

    TFile* fin = TFile::Open(argv[1]);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input file." << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)fin->Get("hepmc3_tree");
    if (!tree) {
        std::cerr << "Tree 'hepmc3_tree' not found." << std::endl;
        return 1;
    }

    HepMC3::GenEventData* sample = nullptr;
    tree->SetBranchAddress("hepmc3_event", &sample);

    const int nEvents = tree->GetEntries();
    analyzer.set_sample_count(nEvents);  
    for (int ie = 0; ie < nEvents; ++ie) {
        tree->GetEntry(ie);
        for (const auto& p : sample->particles) {
            smash::PdgCode pdg(std::to_string(p.pid));
            if (pdg == smash::PdgCode::invalid()) continue;

            double pt = p.momentum.pt();
            double eta = p.momentum.eta();
            double phi = p.momentum.phi();
            double y = p.momentum.rap();
            if (phi > M_PI) phi -= 2 * M_PI;
            if (phi < -M_PI) phi += 2 * M_PI;
            bool is_charged = std::abs(pdg.charge()) > 1e-4;
            if (is_charged){
            total += 1.0;
            if (std::abs(eta) < 0.5) {
                dn_deta += 1.0;
                mean_pt += pt;}
            }
            analyzer.fill(p.pid, eta, pt, phi, y, is_charged);
        }
    }
    std::cout << "Total charged particles: " << total/200. << std::endl;
    std::cout << "dN/deta: " << dn_deta/200. << std::endl;
    std::cout << "Mean pT: " << mean_pt/dn_deta << std::endl;
    TFile* fout = TFile::Open(argv[2], "RECREATE");
    analyzer.write(fout);
    fout->Close();
    fin->Close();

    return 0;
}
