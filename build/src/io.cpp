#include "io.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

void rappas::io::create_directory(const std::string& path)
{
    /// Create if does not exist
    if (!fs::is_directory(path) || !fs::exists(path))
    {
        if (!fs::create_directory(path))
        {
            throw std::runtime_error("Cannot create directory " + path);
        }
    }
}
