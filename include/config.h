#ifndef CONFIG_H
#define CONFIG_H


#include <vector>
#include <string>

struct Config {
    int n_max;
    std::vector<int> pids;
    int pt_bins = 60;
    double max_pt = 3.0;
    int eta_bins = 50;
    double max_eta = 5.0;
    std::vector<double> pt_bins_custom;
    std::vector<double> eta_bins_custom;
    bool has_custom_pt_bins = false;
    bool has_custom_eta_bins = false;
    bool calculate_charged = true;
    std::string input_type = "smash_hepmc3";
    std::string output_type = "root";  // default fallback

    const std::vector<double> pt_grid_values = {
        0.0058285707, 0.0307765083, 0.0759334473, 0.1417947871, 0.2291019148,
        0.3388993645, 0.4726129229, 0.6321707327, 0.8201969196, 1.0403351122,
        1.2978246292, 1.6006235312, 1.9619037526, 2.4068597001, 3.0
    };
    
    const std::vector<double> pt_grid_weights = {
        0.0149655027, 0.0349872889, 0.0554079050, 0.0764374852, 0.0983495971,
        0.1214809798, 0.1462643912, 0.1732834496, 0.2033671894, 0.2377652687,
        0.2785035070, 0.3291972095, 0.3972568279, 0.5017140708, 0.7200944214
    };
    
    const std::vector<double> phi_grid_values = {
        0.0215871423, 0.1131855286, 0.2757236754, 0.5054289394, 0.7969218191,
        1.1433710829, 1.5366566333, 1.9675603507, 2.4259822971, 2.9011774970,
        3.3820078102, 3.8572030101, 4.3156249565, 4.7465286739, 5.1398142242,
        5.4862634881, 5.7777563678, 6.0074616317, 6.1699997786, 6.2615981649
    };
    
    const std::vector<double> phi_grid_weights = {
        0.0553360354, 0.1275531536, 0.1968900466, 0.2616215996, 0.3202229156,
        0.3713190733, 0.4137120591, 0.4464080931, 0.4686407584, 0.4798889188,
        0.4798889188, 0.4686407584, 0.4464080931, 0.4137120591, 0.3713190733,
        0.3202229156, 0.2616215996, 0.1968900466, 0.1275531536, 0.0553360354
    };


};



Config read_config(const std::string& filename);

#endif // CONFIG_H