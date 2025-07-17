#include "config.h"


Config read_config(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename);
    Config cfg;

    cfg.n_max = config["n_max"].as<int>();
    cfg.pids = config["pids"].as<std::vector<int>>();

    cfg.pt_bins = config["pt_bins"] ? config["pt_bins"].as<int>() : 60;
    cfg.max_pt = config["max_pt"] ? config["max_pt"].as<double>() : 3.0;
    cfg.eta_bins = config["eta_bins"] ? config["eta_bins"].as<int>() : 50;
    cfg.max_eta = config["max_eta"] ? config["max_eta"].as<double>() : 5.0;

    if (config["pt_bins_custom"]) {
        for (const auto& v : config["pt_bins_custom"])
            cfg.pt_bins_custom.push_back(v.as<double>());
        cfg.has_custom_pt_bins = true;
    }

    if (config["eta_bins_custom"]) {
        for (const auto& v : config["eta_bins_custom"])
            cfg.eta_bins_custom.push_back(v.as<double>());
        cfg.has_custom_eta_bins = true;
    }

    cfg.calculate_charged = config["calculate_charged"] ?
                            config["calculate_charged"].as<bool>() : true;

    cfg.input_type = config["input_type"] ?
                     config["input_type"].as<std::string>() : "smash_hepmc3";

    cfg.output_type = config["output_type"] ?
                      config["output_type"].as<std::string>(): "root";  //

    return cfg;
}
