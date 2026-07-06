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

        // Build set of incoming particle 1-based indices (particles added as
        // incoming to the main IP vertex, i.e. links2[j] == -1). In List modus
        // these are the initial input particles which share status==1 with the
        // final state particles and must be excluded explicitly.
        std::unordered_set<int> incoming_indices;
        for (size_t i = 0; i < sample->links1.size(); ++i) {
            if (sample->links2[i] == -1 && sample->links1[i] > 0) {
                incoming_indices.insert(sample->links1[i]);
            }
        }

        for (size_t i = 0; i < sample->particles.size(); ++i) {
            const auto& p = sample->particles[i];

            // Skip incoming (initial state) particles
            if (incoming_indices.count((int)i + 1)) continue;

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
        if (pdg == smash::PdgCode::invalid()){
            std::cout << " invalid particle found" << std::endl;
            continue;
        }

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
        //print if it is charged ;
        std::cout << "Processing PID " << pid << " (charged: " << is_charged << ", requested: " << is_requested << ")" << std::endl;

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



void read_oscar_sampler(const std::string& filename, qvector_writter& analyzer,
                        double& total, double& dn_deta, double& mean_pt, int& nEvents) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Cannot open OSCAR sampler file: " << filename << std::endl;
        std::exit(1);
    }

    std::string line;
    bool in_event = false;
    nEvents = 0;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        // Header lines
        if (line[0] == '#') {
            // "# event <i>"
            if (line.find("# event") != std::string::npos) in_event = true;
            continue;
        }

        // End marker: " end <i>"
        {
            std::istringstream iss(line);
            std::string first;
            if (iss >> first) {
                if (first == "end") {
                    in_event = false;
                    ++nEvents;
                    continue;
                }
            }
        }

        if (!in_event) continue;

        // Particle line:
        // pid t x y z mass E px py pz
        std::istringstream iss(line);
        int pid;
        double t, x, y, z, mass, E, px, py, pz;

        if (!(iss >> pid >> t >> x >> y >> z >> mass >> E >> px >> py >> pz)) {
            // ignore malformed lines
            continue;
        }

        smash::PdgCode pdg(std::to_string(pid));
        if (pdg == smash::PdgCode::invalid()) continue;

        double pt  = std::sqrt(px * px + py * py);
        double p   = std::sqrt(pt * pt + pz * pz);

        double y_rap = 0.5 * std::log((E + pz) / (E - pz + 1e-10));
        double eta   = 0.5 * std::log((p + pz) / (p - pz + 1e-10));
        double phi   = std::atan2(py, px);

        if (phi > M_PI)  phi -= 2.0 * M_PI;
        if (phi < -M_PI) phi += 2.0 * M_PI;

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

void read_iss_oscar(const std::string& filename, qvector_writter& analyzer,
                    double& total, double& dn_deta, double& mean_pt, int& nEvents)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Cannot open iSS OSCAR file: " << filename << std::endl;
        std::exit(1);
    }

    // Reset
    nEvents  = 0;

    // --- Skip the 2 run-header lines (can be more in some generators, but for your file it's 2)
    // Example:
    //   final_id_p_x
    //   3DHydro 1.1 (197,79)+(197,79) ...
    std::string line;
    if (!std::getline(fin, line)) { fin.close(); return; }
    if (!std::getline(fin, line)) { fin.close(); return; }

    // Now parse repeated blocks:
    //   event_id  N  0  0
    //   (N lines): i pid px py pz E m x y z t
    while (true) {
        int event_id = 0;
        int N = 0;
        int dummy1 = 0, dummy2 = 0;

        // Read event header (skip blank lines if any)
        bool got_header = false;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            if (iss >> event_id >> N >> dummy1 >> dummy2) {
                got_header = true;
                break;
            }
            // If we can't parse it, just keep scanning (some files have extra text)
        }

        if (!got_header) break;  // EOF

        ++nEvents;

        for (int ip = 0; ip < N; ++ip) {
            if (!std::getline(fin, line)) break;
            if (line.empty()) { --ip; continue; }  // be robust to blank lines

            std::istringstream iss(line);

            int idx = 0;
            int pid = 0;
            double px = 0.0, py = 0.0, pz = 0.0, E = 0.0;
            double m  = 0.0;
            double x  = 0.0, yx = 0.0, z = 0.0, t = 0.0;

            // Format:
            // idx pid px py pz E m x y z t
            if (!(iss >> idx >> pid >> px >> py >> pz >> E >> m >> x >> yx >> z >> t)) {
                // malformed line; skip it safely
                continue;
            }

            smash::PdgCode pdg(std::to_string(pid));
            if (pdg == smash::PdgCode::invalid()) {
                // keep going, but warn once in a while if you want
                continue;
            }

            const double pt  = std::sqrt(px * px + py * py);
            const double p   = std::sqrt(pt * pt + pz * pz);

            // robust rapidity/eta (avoid division by zero)
            const double y_rap = 0.5 * std::log((E + pz) / (E - pz + 1e-12));
            const double eta   = 0.5 * std::log((p + pz) / (p - pz + 1e-12));
            double phi         = std::atan2(py, px);

            if (phi > M_PI)  phi -= 2.0 * M_PI;
            if (phi < -M_PI) phi += 2.0 * M_PI;

            const bool is_charged = std::abs(pdg.charge()) > 1e-4;

            // Your existing convention: count charged particles, and use |eta|<0.5 for dn/deta and mean_pt
            if (is_charged) {
                total += 1.0;
                if (std::abs(eta) < 0.5) {
                    dn_deta += 1.0;
                    mean_pt += pt;
                }
            }

            analyzer.fill(pid, eta, pt, phi, y_rap, is_charged);
        }
    }

    fin.close();
    analyzer.set_sample_count(nEvents);
}


