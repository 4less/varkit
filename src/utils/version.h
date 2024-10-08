//
// Created by fritsche on 05/07/22.
//

#pragma once

#include <string>
#include <tuple>

//using namespace constants;

// This uint32_t for version may not be changed, otherwise some data
// will not be consistent anymore and crash.
// It is vital that the versioning class in its content and struct size
// stays untouched.
using VersionType = uint32_t;
using VersionMajorType = VersionType;
using VersionMinorType = VersionType;
using VersionPatchType = VersionType;
using VersionTweakType = VersionType;

class Version {
    const VersionMajorType version_major = 0;
    const VersionMinorType version_minor = 0;
    const VersionPatchType version_patch = 0;
    const VersionTweakType version_tweak = 0;
    using VersionTuple = std::tuple<VersionMajorType, VersionMinorType, VersionPatchType,  VersionTweakType>;



    static VersionTuple FromString(std::string version_str) noexcept {
        std::string token = "";
        VersionMajorType version_major = 0;
        VersionMinorType version_minor = 0;
        VersionPatchType version_patch = 0;
        VersionTweakType version_tweak = 0;
        int state = 0;
        bool valid = false;

        auto set_version = [&](int type, size_t value) {
            if (type == 0) {
                version_major = value;
            } else if (type == 1) {
                version_minor = value;
            } else if (type == 2) {
                version_patch = value;
            } else if (type == 3) {
                version_tweak = value;
            }
        };

        for (auto i = 0; i < version_str.length() + 1; i++) {
            if (valid && (i == version_str.length() || version_str[i] == '.')) {
                if (state < 4) {
                    set_version(state, std::stoul(token));
                    token = "";
                    state++;
                } else {
                    return { 0, 0, 0, 0 };
                }
            } else if (std::isdigit(version_str[i])) { valid = true; token += version_str[i]; }
            else {
                return { 0, 0, 0, 0 };
            }
        }
        return { version_major, version_minor, version_patch, version_tweak };
    }

public:
    Version() {};

    Version(std::string version_str) :
            version_major(std::get<0>(FromString(version_str))),
            version_minor(std::get<1>(FromString(version_str))),
            version_patch(std::get<2>(FromString(version_str))),
            version_tweak(std::get<3>(FromString(version_str))) {};

    Version(VersionMajorType major, VersionMinorType minor, VersionPatchType patch, VersionTweakType tweak) :
            version_major(major),
            version_minor(minor),
            version_patch(patch),
            version_tweak(tweak) {};

    Version(VersionMajorType major, VersionMinorType minor, VersionPatchType patch) :
            version_major(major),
            version_minor(minor),
            version_patch(patch) {};

    Version(VersionMajorType major, VersionMinorType minor) :
            version_major(major),
            version_minor(minor) {};

    Version(VersionMajorType major) :
            version_major(major) {};

    bool operator==(const Version &other_version) const {
        return
                version_major == other_version.version_major &&
                version_minor == other_version.version_minor &&
                version_patch == other_version.version_patch &&
                version_tweak == other_version.version_tweak;
    }
    bool operator<(const Version &other_version) const
    {
        if (version_major != other_version.version_major)
            return version_major < other_version.version_major;
        if (version_minor != other_version.version_minor)
            return version_minor < other_version.version_minor;
        if (version_patch != other_version.version_patch)
            return version_patch < other_version.version_patch;
        return version_tweak < other_version.version_tweak;
    }
    bool operator>(const Version &other_version) const
    {
        return other_version < *this;
    }
    bool operator<=(const Version &other_version) const {
        return *this < other_version || *this == other_version;
    }
    bool operator>=(const Version &other_version) const {
        return *this > other_version || *this == other_version;
    }

    std::string ToString() const {
        return
            std::to_string(version_major) + "." +
            std::to_string(version_minor) + "." +
            std::to_string(version_patch) + "." +
            std::to_string(version_tweak);
    }
};