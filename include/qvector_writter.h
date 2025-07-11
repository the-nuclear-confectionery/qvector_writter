#ifndef QVECTOR_WRITTER_H
#define QVECTOR_WRITTER_H

#include <map>
#include <string>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

// ==== Configuration Struct ====
struct Config {
    int n_max;
    std::vector<int> pids;  // ✅ Restored — for identified particles

    int eta_bins;
    double max_eta;
    int pt_bins;
    double max_pt;

    bool calculate_charged;

    bool has_custom_eta_bins = false;
    bool has_custom_pt_bins = false;
    std::vector<double> eta_bins_custom;
    std::vector<double> pt_bins_custom;
};

// ==== Read config from YAML ====
Config read_config(const std::string& filename);

// ==== Q-vector Writer Class ====
class qvector_writter {
public:
    qvector_writter(const Config& config);

    void fill(int pid, double eta, double pt, double phi, double y, bool is_charged);
    void write(TFile* outfile);

    // ✅ New function to set the sample count
    void set_sample_count(int nsamples);

private:
    Config config_;

    // Histograms for identified particles
    std::map<int, std::vector<TH2D>> hReQ_;
    std::map<int, std::vector<TH2D>> hImQ_;

    // Histograms for charged particles
    std::vector<TH2D> hReQ_charged_;
    std::vector<TH2D> hImQ_charged_;

    // ✅ Histogram to store the number of samples
    TH1D hSampleCounter_;
};

#endif
