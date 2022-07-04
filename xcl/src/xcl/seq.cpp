#include "xcl/seq.h"

using namespace xcl;

#ifdef SEQ_TYPE_DNA
const seq_traits_impl<dna>::char_type seq_traits_impl<dna>::code_to_key[] = {'A', 'C', 'G', 'T'};

#elif SEQ_TYPE_AA

const seq_traits_impl<aa>::char_type seq_traits_impl<aa>::code_to_key[] = {
            'R', 'H', 'K',
            'D', 'E',
            'S', 'T', 'N', 'Q',
            'C', 'G', 'P',
            'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V'
        };
#endif
