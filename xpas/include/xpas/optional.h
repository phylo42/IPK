#ifndef XPAS_OPTIONAL_H
#define XPAS_OPTIONAL_H

#if __has_include(<optional>)
#include <optional>
using std::optional;
using std::nullopt;
#elif __has_include(<experimental/optional>)
#include <experimental/optional>
using std::experimental::optional;
using std::experimental::nullopt;
#endif

#endif //XPAS_OPTIONAL_H
