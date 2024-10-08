//
// Created by fritsche on 22/03/2022.
//

#pragma once

#include <OptionsContainer.h>

namespace ProfileMode {
    static void Profile(AbundanceEstimator &ae, VarkitOptionsContainer &options, Sample& sample) {
        ae.Clear();
        if (options.raw_taxa) {
            ae.SetRawOutput(options.RawTaxaOutputFile(sample));
        }
        auto classification_output = options.ClassificationOutputFile(sample);
        if (!Utils::exists(classification_output)) {
            std::cerr << classification_output << " does not exist" << std::endl;
            exit(9);
        }
        auto profile_output = options.ProfileOutputFile(sample);
        std::cout << "Load read classifications" << std::endl;
        ae.ReadClassificationFile(classification_output.c_str());
//        ae.DebugFunction();

        std::cout << "Profile" << std::endl;
        ae.Profile(options.ProfileOutputFile(sample));
    }

    static void Profile(AbundanceEstimator &ae, VarkitOptionsContainer &options, std::string classification_file) {
        ae.Clear();
        ae.ReadClassificationFile(classification_file.c_str());
        ae.Profile(options.ProfileOutputFileFromClass(classification_file));
    }

    static void Run(AbundanceEstimator &ae, VarkitOptionsContainer &options, Sample& sample) {
        Profile(ae, options, sample);
    }
}