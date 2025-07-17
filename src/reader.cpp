#include "reader.h"


void read_smash_hepmc3(const std::string& filename, qvector_writter& analyzer,
                     double& total, double& dn_deta, double& mean_pt, int& nEvents) {
    TFile* fin = TFile::Open(filename.c_str());
    if (!fin || fin->IsZombie()) {
        std::cerr << "Cannot open input ROOT file." << std::endl;
        exit(1);
    }

    TTree* tree = (TTree*)fin->Get("hepmc3_tree");
    if (!tree) {
        std::cerr << "Tree 'hepmc3_tree' not found." << std::endl;
        exit(1);
    }

    HepMC3::GenEventData* sample = nullptr;
    tree->SetBranchAddress("hepmc3_event", &sample);

    nEvents = tree->GetEntries();
    analyzer.set_sample_count(nEvents);

    for (int ie = 0; ie < nEvents; ++ie) {
        tree->GetEntry(ie);

        // Determine the last incoming particle index
        int max_incoming_particle_index = -1;
        for (size_t i = 0; i < sample->links1.size(); ++i) {
            if (sample->links2[i] == -1 && sample->links1[i] > max_incoming_particle_index) {
                max_incoming_particle_index = sample->links1[i];
            }
        }

        for (size_t i = 0; i < sample->particles.size(); ++i) {
            const auto& p = sample->particles[i];

            // Skip incoming particles
            if ((int)i < max_incoming_particle_index) continue;

            // Only final state particles
            if (p.status != 1) continue;

            smash::PdgCode pdg(std::to_string(p.pid));
            if (pdg == smash::PdgCode::invalid()) continue;

            double pt = p.momentum.pt();
            double eta = p.momentum.eta();
            double phi = p.momentum.phi();
            double y = p.momentum.rap();

            if (phi > M_PI) phi -= 2 * M_PI;
            if (phi < -M_PI) phi += 2 * M_PI;

            bool is_charged = std::abs(pdg.charge()) > 1e-4;
            if (is_charged) {
                total += 1.0;
                if (std::abs(eta) < 0.5) {
                    dn_deta += 1.0;
                    mean_pt += pt;
                }
            }

            analyzer.fill(p.pid, eta, pt, phi, y, is_charged);
        }
    }

    fin->Close();
}

void read_oscar_file(const std::string& filename, qvector_writter& analyzer,
                     double& total, double& dn_deta, double& mean_pt, int& nEvents) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Cannot open OSCAR file: " << filename << std::endl;
        exit(1);
    }

    std::string line;
    bool in_event = false;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') {
            if (line.find("# event") != std::string::npos && line.find("out") != std::string::npos)
                in_event = true;
            else if (line.find("# event") != std::string::npos && line.find("end") != std::string::npos) {
                in_event = false;
                ++nEvents;
            }
            continue;
        }

        if (!in_event) continue;

        std::istringstream iss(line);
        double t, x, y, z, mass, p0, px, py, pz;
        int pid, local_id, charge_int;

        if (!(iss >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pid >> local_id >> charge_int))
            continue;

        double pt = std::sqrt(px * px + py * py);
        double p = std::sqrt(pt * pt + pz * pz);
        double E = p0;
        double y_rap = 0.5 * std::log((E + pz) / (E - pz + 1e-10));
        double eta = 0.5 * std::log((p + pz) / (p - pz + 1e-10));
        double phi = std::atan2(py, px);

        if (phi > M_PI) phi -= 2 * M_PI;
        if (phi < -M_PI) phi += 2 * M_PI;


        smash::PdgCode pdg(std::to_string(pid));
        if (pdg == smash::PdgCode::invalid()) continue;

        bool is_charged = std::abs(pdg.charge()) > 1e-4;
        if (is_charged) {
            total += 1.0;
            if (std::abs(eta) < 0.5) {
                dn_deta += 1.0;
                mean_pt += pt;
            }
        }

        analyzer.fill(pid, eta, pt, phi, y_rap, is_charged);
    }

    fin.close();
    analyzer.set_sample_count(nEvents);
}



void read_afterdecays(const std::string& input_filename,
                      const Config& cfg,
                      qvector_writter& analyzer,
                      double& total, double& dn_deta, double& mean_pt, int& nEvents) {
    std::ifstream input(input_filename);
    if (!input.is_open()) {
        std::cerr << "Cannot open afterdecays .dat file: " << input_filename << std::endl;
        exit(1);
    }

    const auto& pt_grid = cfg.pt_grid_values;
    const auto& phi_grid = cfg.phi_grid_values;
    const auto& phi_weights = cfg.phi_grid_weights;

    const size_t nPT = pt_grid.size();
    const size_t nPhi = phi_grid.size();

    std::map<int, std::vector<double>> spectra_data_flat;
    std::string line;

    while (std::getline(input, line)) {
        if (line.empty()) continue;

        int pid = std::stoi(line);
        std::vector<double> grid_flat(nPT * nPhi, 0.0);

        for (size_t iPhi = 0; iPhi < nPhi; ++iPhi) {
            if (!std::getline(input, line)) break;
            std::istringstream iss(line);
            for (size_t iPT = 0; iPT < nPT; ++iPT) {
                double val;
                if (!(iss >> val)) {
                    std::cerr << "Incomplete row for PID " << pid << std::endl;
                    exit(1);
                }
                grid_flat[iPT * nPhi + iPhi] = val;  // Transpose: [pt][phi]
            }
        }

        spectra_data_flat[pid] = std::move(grid_flat);
    }

    input.close();

    // Analyze PID presence
    std::vector<int> used_pids;
    std::vector<int> missing_pids;
    std::vector<int> not_used_pids;

    for (const int pid : cfg.pids) {
        if (spectra_data_flat.find(pid) != spectra_data_flat.end()) {
            used_pids.push_back(pid);
        } else {
            missing_pids.push_back(pid);
        }
    }

    for (const auto& [pid, _] : spectra_data_flat) {
        if (std::find(cfg.pids.begin(), cfg.pids.end(), pid) == cfg.pids.end()) {
            not_used_pids.push_back(pid);
        }
    }

    if (used_pids.empty()) {
        std::cerr << "Warning: None of the requested PIDs were found in file: " << input_filename << std::endl;
        return;
    }

    if (!missing_pids.empty()) {
        std::cerr << "Warning: Missing PIDs in file: ";
        for (int pid : missing_pids) std::cerr << pid << " ";
        std::cerr << std::endl;
    }

    if (!not_used_pids.empty()) {
        std::cerr << "Info: Unused PIDs in file: ";
        for (int pid : not_used_pids) std::cerr << pid << " ";
        std::cerr << std::endl;
    }

    // Process all PIDs
    for (const auto& [pid, grid_flat] : spectra_data_flat) {
        smash::PdgCode pdg(std::to_string(pid));
        bool is_charged = std::abs(pdg.charge()) > 1e-4;
        bool is_requested = std::find(cfg.pids.begin(), cfg.pids.end(), pid) != cfg.pids.end();

        // Reconstruct 2D grid: [phi][pt]
        std::vector<std::vector<double>> grid_2D(nPhi, std::vector<double>(nPT, 0.0));
        for (size_t iPT = 0; iPT < nPT; ++iPT) {
            for (size_t iPhi = 0; iPhi < nPhi; ++iPhi) {
                grid_2D[iPhi][iPT] = grid_flat[iPT * nPhi + iPhi];
            }
        }

        // Integrate for total/dn_deta/mean_pt
        if (is_charged) {
            for (size_t iPT = 0; iPT < nPT; ++iPT) {
                double pt = pt_grid[iPT];
                for (size_t iPhi = 0; iPhi < nPhi; ++iPhi) {
                    double val = pt * grid_2D[iPhi][iPT];
                    if (val > 0.0) {
                        total += val;
                        dn_deta += val;
                        mean_pt += pt * val;
                    }
                }
            }
        }

        if (is_charged || is_requested) {
            analyzer.fill_afterdecays(pid, pt_grid, phi_grid, phi_weights, grid_2D, is_charged);
        }
    }

    ++nEvents;
    analyzer.set_sample_count(nEvents);
}


