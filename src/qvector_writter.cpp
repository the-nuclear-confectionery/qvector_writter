#include "qvector_writter.h"
#include "TH1D.h"
#include "TString.h"

#include <cmath>
#include <iostream>

// ==== Constructor ====
qvector_writter::qvector_writter(const Config& config) :
    config_(config),
    hSampleCounter_("hSampleCounter", "Number of samples", 1, 0, 1) {

    for (int pid : config_.pids) {
        std::vector<TH2D> hre, him;
        for (int n = 0; n <= config_.n_max; ++n) {
            TString name_re = TString::Format("ReQ_pid%d_n%d", pid, n);
            TString name_im = TString::Format("ImQ_pid%d_n%d", pid, n);
            TString title   = TString::Format("Q_{%d} for pid %d", n, pid);

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

    if (config_.calculate_charged) {
        for (int n = 0; n <= config_.n_max; ++n) {
            TString name_re = TString::Format("ReQ_charged_n%d", n);
            TString name_im = TString::Format("ImQ_charged_n%d", n);
            TString title   = TString::Format("Q_{%d} for charged", n);

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

// ==== Set sample count ====
void qvector_writter::set_sample_count(int nsamples) {
    hSampleCounter_.Fill(0.5, nsamples);
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

void qvector_writter::fill_afterdecays(
    int pid,
    const std::vector<double>& pt_grid,                 // size: nPT
    const std::vector<double>& phi_grid,                // size: nPhi
    const std::vector<double>& phi_weights,             // size: nPhi
    const std::vector<std::vector<double>>& grid,       // [nPhi][nPT]
    bool is_charged)
{
    const int n_max = config_.n_max;
    const size_t nPT = pt_grid.size();
    const size_t nPhi = phi_grid.size();

    // Step 1: Compute dQn/dpT at each pt_grid point
    std::vector<std::vector<double>> Qn_real(n_max + 1, std::vector<double>(nPT, 0.0));
    std::vector<std::vector<double>> Qn_imag(n_max + 1, std::vector<double>(nPT, 0.0));

    for (size_t iPT = 0; iPT < nPT; ++iPT) {
        double pt = pt_grid[iPT];
        for (int n = 0; n <= n_max; ++n) {
            double Re = 0.0, Im = 0.0;
            for (size_t iPhi = 0; iPhi < nPhi; ++iPhi) {
                double phi = phi_grid[iPhi];
                double weight = phi_weights[iPhi];
                double val = grid[iPhi][iPT];  // dN / (pT dpT dphi)
                Re += weight * pt * val * std::cos(n * phi);
                Im += weight * pt * val * std::sin(n * phi);
            }
            Qn_real[n][iPT] = Re;
            Qn_imag[n][iPT] = Im;
        }
    }

    // Step 2: Interpolate dQn/dpT onto histogram bin centers and multiply by dpT to get Qn(pt)
    for (int n = 0; n <= n_max; ++n) {
        // Create PCHIP interpolators
        auto spline_real = boost::math::interpolators::pchip(
            std::vector<double>(pt_grid.begin(), pt_grid.end()),
            std::vector<double>(Qn_real[n].begin(), Qn_real[n].end()));

        auto spline_imag = boost::math::interpolators::pchip(
            std::vector<double>(pt_grid.begin(), pt_grid.end()),
            std::vector<double>(Qn_imag[n].begin(), Qn_imag[n].end()));

        for (int ipt = 0; ipt < hReQ_charged_[n].GetNbinsY(); ++ipt) {
            double pt = hReQ_charged_[n].GetYaxis()->GetBinCenter(ipt + 1);
            double dpT = hReQ_charged_[n].GetYaxis()->GetBinWidth(ipt + 1);

            double Re_interp = 0.0;
            double Im_interp = 0.0;

            // Only interpolate within valid pt range
            if (pt >= pt_grid.front() && pt <= pt_grid.back()) {
                Re_interp = spline_real(pt);
                Im_interp = spline_imag(pt);

                // Convert from dQ/dpT to Q by multiplying by dpT
                Re_interp *= dpT;
                Im_interp *= dpT;
            }

            if (config_.calculate_charged && is_charged) {
                hReQ_charged_[n].Fill(0.0, pt, Re_interp);
                hImQ_charged_[n].Fill(0.0, pt, Im_interp);
            }

            if (hReQ_.count(pid) > 0) {
                hReQ_[pid][n].Fill(0.0, pt, Re_interp);
                hImQ_[pid][n].Fill(0.0, pt, Im_interp);
            }
        }
    }
}

// ==== Write to ROOT file ====
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

    hSampleCounter_.Write();
}

void qvector_writter::write_text(std::ostream& out) const {
    auto write_hist = [&](const TH2D& h) {
        out << "# Histogram: " << h.GetName() << "\n";
        out << "# X (eta/y)\tY (pt)\tContent\n";
        for (int ix = 1; ix <= h.GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= h.GetNbinsY(); ++iy) {
                double x = h.GetXaxis()->GetBinCenter(ix);
                double y = h.GetYaxis()->GetBinCenter(iy);
                double z = h.GetBinContent(ix, iy);
                out << x << "\t" << y << "\t" << z << "\n";
            }
        }
        out << "\n";
    };

    for (const auto& [pid, vec] : hReQ_)
        for (const auto& h : vec) write_hist(h);
    for (const auto& [pid, vec] : hImQ_)
        for (const auto& h : vec) write_hist(h);

    for (const auto& h : hReQ_charged_) write_hist(h);
    for (const auto& h : hImQ_charged_) write_hist(h);
}

