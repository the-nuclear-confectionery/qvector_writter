#include "qvector_writter.h"
#include <TFile.h>
#include <TString.h>
#include <cmath>
#include <yaml-cpp/yaml.h>
#include <iostream>

qvector_writter::qvector_writter(const Config& config) : config_(config),
    hSampleCounter_("hSampleCounter", "Number of samples", 1, 0, 1)  // ✅ Initialize sample counter
{
    // Identified particles
    for (auto pid : config_.pids) {
        std::vector<TH2D> hre, him;
        for (int n = 0; n <= config_.n_max; ++n) {
            TString name_re = TString::Format("ReQ_pid%d_n%d", pid, n);
            TString name_im = TString::Format("ImQ_pid%d_n%d", pid, n);
            TString title = TString::Format("Q_{%d} for pid %d", n, pid);

            if (config_.has_custom_eta_bins && config_.has_custom_pt_bins) {
                hre.emplace_back(name_re, title,
                    config_.eta_bins_custom.size() - 1, config_.eta_bins_custom.data(),
                    config_.pt_bins_custom.size() - 1, config_.pt_bins_custom.data());
                him.emplace_back(name_im, title,
                    config_.eta_bins_custom.size() - 1, config_.eta_bins_custom.data(),
                    config_.pt_bins_custom.size() - 1, config_.pt_bins_custom.data());
            } else {
                hre.emplace_back(name_re, title,
                    config_.eta_bins, -config_.max_eta, config_.max_eta,
                    config_.pt_bins, 0, config_.max_pt);
                him.emplace_back(name_im, title,
                    config_.eta_bins, -config_.max_eta, config_.max_eta,
                    config_.pt_bins, 0, config_.max_pt);
            }
        }
        hReQ_[pid] = hre;
        hImQ_[pid] = him;
    }

    // Charged particles
    if (config_.calculate_charged) {
        for (int n = 0; n <= config_.n_max; ++n) {
            TString name_re = TString::Format("ReQ_charged_n%d", n);
            TString name_im = TString::Format("ImQ_charged_n%d", n);
            TString title = TString::Format("Q_{%d} for charged", n);

            if (config_.has_custom_eta_bins && config_.has_custom_pt_bins) {
                hReQ_charged_.emplace_back(name_re, title,
                    config_.eta_bins_custom.size() - 1, config_.eta_bins_custom.data(),
                    config_.pt_bins_custom.size() - 1, config_.pt_bins_custom.data());
                hImQ_charged_.emplace_back(name_im, title,
                    config_.eta_bins_custom.size() - 1, config_.eta_bins_custom.data(),
                    config_.pt_bins_custom.size() - 1, config_.pt_bins_custom.data());
            } else {
                hReQ_charged_.emplace_back(name_re, title,
                    config_.eta_bins, -config_.max_eta, config_.max_eta,
                    config_.pt_bins, 0, config_.max_pt);
                hImQ_charged_.emplace_back(name_im, title,
                    config_.eta_bins, -config_.max_eta, config_.max_eta,
                    config_.pt_bins, 0, config_.max_pt);
            }
        }
    }
}

void qvector_writter::set_sample_count(int nsamples) {
    hSampleCounter_.Fill(0.5, nsamples);  // ✅ Added
}

void qvector_writter::fill(int pid, double eta, double pt, double phi, double y, bool is_charged) {
    for (int n = 0; n <= config_.n_max; ++n) {
        double cos_nphi = (n == 0) ? 1.0 : std::cos(n * phi);
        double sin_nphi = (n == 0) ? 0.0 : std::sin(n * phi);

        // Charged
        if (config_.calculate_charged && is_charged) {
            hReQ_charged_[n].Fill(y, pt, cos_nphi);
            hImQ_charged_[n].Fill(y, pt, sin_nphi);
        }

        // Identified
        if (hReQ_.count(pid) > 0) {
            hReQ_[pid][n].Fill(eta, pt, cos_nphi);
            hImQ_[pid][n].Fill(eta, pt, sin_nphi);
        }
    }
}

void qvector_writter::write(TFile* outfile) {
    outfile->cd();

    for (auto& [pid, vec] : hReQ_)
        for (auto& h : vec) h.Write();
    for (auto& [pid, vec] : hImQ_)
        for (auto& h : vec) h.Write();

    for (auto& h : hReQ_charged_)
        h.Write();
    for (auto& h : hImQ_charged_)
        h.Write();

    hSampleCounter_.Write();  // ✅ Added
}

Config read_config(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename);
    Config cfg;
    cfg.n_max = config["n_max"].as<int>();
    cfg.pids = config["pids"].as<std::vector<int>>();

    cfg.pt_bins = config["pt_bins"].as<int>(60);
    cfg.max_pt = config["max_pt"].as<double>(3.0);
    cfg.eta_bins = config["eta_bins"].as<int>(50);
    cfg.max_eta = config["max_eta"].as<double>(5.0);

    if (config["pt_bins_custom"]) {
        for (const auto& v : config["pt_bins_custom"])
            cfg.pt_bins_custom.push_back(v.as<double>());
    }

    if (config["eta_bins_custom"]) {
        for (const auto& v : config["eta_bins_custom"])
            cfg.eta_bins_custom.push_back(v.as<double>());
    }

    return cfg;
}
