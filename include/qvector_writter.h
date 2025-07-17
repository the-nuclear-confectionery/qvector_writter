#ifndef QVECTOR_WRITTER_H
#define QVECTOR_WRITTER_H

#include "TH2D.h"
#include "TFile.h"
#include <map>
#include <vector>
#include <string>
#include "config.h"
#include <boost/math/interpolators/pchip.hpp>

#include <complex>

// Forward declare Config
struct Config;

class qvector_writter {
public:
    qvector_writter(const Config& config);

    void fill(int pid, double eta, double pt, double phi, double y, bool is_charged);

    void fill_afterdecays(  int pid,
    const std::vector<double>& pt_grid,                 // size: nPT
    const std::vector<double>& phi_grid,                // size: nPhi
    const std::vector<double>& phi_weights,             // size: nPhi
    const std::vector<std::vector<double>>& grid,       // [nPhi][nPT]
    bool is_charged);

    void set_sample_count(int nsamples);
    void write(TFile* outfile);
    void write_text(std::ostream& out) const;


private:
    const Config& config_;
    std::map<int, std::vector<TH2D>> hReQ_, hImQ_;
    std::vector<TH2D> hReQ_charged_, hImQ_charged_;
    TH1D hSampleCounter_;
};


#endif // QVECTOR_WRITTER_H