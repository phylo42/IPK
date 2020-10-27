#ifndef SHERPAS_OPTIONAL_H
#define SHERPAS_OPTIONAL_H

#ifdef __APPLE__

/// Seriously Apple? *Sigh*
#include <experimental/optional>
using std::experimental::optional;
using std::experimental::nullopt;

#else

#include <optional>
using std::optional;
using std::nullopt;

#endif

#endif //SHERPAS_OPTIONAL_H
