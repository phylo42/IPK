#ifndef XPAS_BUILD_SEQ_RECORD_H
#define XPAS_BUILD_SEQ_RECORD_H

#include <string>

namespace xpas
{
    class seq_record
    {
    public:
        seq_record() = delete;
        seq_record(std::string header, std::string sequence) noexcept;
        seq_record(seq_record&& other) noexcept = default;
        seq_record(const seq_record&) = default;
        ~seq_record() = default;

        [[nodiscard]]
        std::string_view header() const noexcept;
        [[nodiscard]]
        std::string_view sequence() const noexcept;

        bool operator==(const seq_record& rhs) const noexcept;
    private:
        std::string _header;
        std::string _sequence;
    };

}
#endif //XPAS_BUILD_SEQ_RECORD_H
