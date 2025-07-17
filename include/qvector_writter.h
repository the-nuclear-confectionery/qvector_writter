#ifndef QVECTOR_WRITTER_H
#define QVECTOR_WRITTER_H

#include "TH2D.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

#include <boost/math/interpolators/pchip.hpp>
#include <map>
#include <vector>
#include <string>
#include "config.h"
#include <cmath>
#include <iostream>
#include <complex>



#include "qvector_writter.h"



// Forward declare Config
struct Config;

class qvector_writter {
public:
    qvector_writter(const Config& config);

    void fill(int pid, double eta, double pt, double phi, double y, bool is_charged);

    void fill_afterdecays(  int pid,
    const std::vector<double>& pt_grid,          
    const std::vector<double>& phi_grid,               
    const std::vector<double>& phi_weights,            
    const std::vector<std::vector<double>>& grid,       
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