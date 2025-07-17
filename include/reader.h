#ifndef READER_H
#define READER_H

#include "config.h"
#include "qvector_writter.h"
#include <string>

void read_smash_hepmc3(const std::string& filename, qvector_writter& analyzer,
                     double& total, double& dn_deta, double& mean_pt, int& nEvents);

void read_oscar_file(const std::string& filename, qvector_writter& analyzer,
                     double& total, double& dn_deta, double& mean_pt, int& nEvents);

void read_afterdecays(const std::string& input_filename,
                      const Config& cfg,
                      qvector_writter& analyzer,
                      double& total, double& dn_deta, double& mean_pt, int& nEvents);


#endif // READER_H
