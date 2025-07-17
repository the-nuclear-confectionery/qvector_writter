#include "config.h"
#include "reader.h"
#include "qvector_writter.h"


int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " input.[root|txt] output.root config.yaml\n";
        return 1;
    }

    const std::string input_file = argv[1];
    const std::string output_file = argv[2];
    const std::string config_file = argv[3];

    // === Load config ===
    Config cfg = read_config(config_file);
    qvector_writter analyzer(cfg);

    double total = 0.;
    double dn_deta = 0.;
    double mean_pt = 0.;
    int nEvents = 0;

    // Read based on input type
    std::cout << "Reading input file: " << input_file << std::endl;
    
    if (cfg.input_type == "smash_hepmc3") {
        read_smash_hepmc3(input_file, analyzer, total, dn_deta, mean_pt, nEvents);
    }
    else if (cfg.input_type == "smash_oscar") {
        read_oscar_file(input_file, analyzer, total, dn_deta, mean_pt, nEvents);
    }
    else if (cfg.input_type == "afterdecays") {
        read_afterdecays(input_file, cfg, analyzer, total, dn_deta, mean_pt, nEvents);
    }
    else {
        std::cerr << "Error: unknown input_type \"" << cfg.input_type << "\" in config file.\n";
        return 1;
    }

    // Output results based on output type
    std::cout << "Writting results to " << output_file << std::endl;

    if (cfg.output_type == "root") {
        std::string root_output_file = output_file + ".root";
        TFile* fout = TFile::Open(root_output_file.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            std::cerr << "Could not open output file: " << root_output_file << std::endl;
            return 1;
        }
        analyzer.write(fout);
        fout->Close();
        std::cout << "Wrote ROOT histograms to: " << root_output_file << std::endl;
    } else if (cfg.output_type == "text") {
        std::string txt_output_file = output_file + ".dat";
        std::ofstream fout_txt(txt_output_file);
        if (!fout_txt.is_open()) {
            std::cerr << "Could not open text output file: " << txt_output_file << std::endl;
            return 1;
        }
        analyzer.write_text(fout_txt);
        fout_txt.close();
        std::cout << "Wrote text histograms to: " << txt_output_file << std::endl;
    } else {
        std::cerr << "Error: unknown output_type \"" << cfg.output_type << "\" in config file.\n";
        return 1;
    }

    return 0;
}
