#include "xpas/seq.h"

using namespace xpas;

#ifdef SEQ_TYPE_DNA
const seq_traits_impl<dna>::char_type seq_traits_impl<dna>::code_to_key[] = {'A', 'C', 'G', 'T'};

#elif SEQ_TYPE_AA

#endif
