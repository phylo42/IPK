#ifndef XPAS_VERSION_H
#define XPAS_VERSION_H

namespace xpas
{
    /// \brief xpas version class. Follows the rules of semantic versioning 2.0.
    /// \sa Semantic versioning 2.0: https://semver.org/
    /// TODO: fill this automatically with the CMakeLists.txt definitions
    struct version
    {
        static constexpr auto major = "0";
        static constexpr auto minor = "1";
        static constexpr auto revision = "7";

        /// \brief Returns the core version number as a string
        static std::string as_string()
        {
            return major + std::string{ "." } + minor + std::string{ "." } + revision;
        }
    };
}

#endif // XPAS_VERSION_H
