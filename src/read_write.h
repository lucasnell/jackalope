#ifndef _READWRITE_H
#define _READWRITE_H


#include <RcppArmadillo.h>

#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "new_variants.h"  // new classes
#include "str_manip.h"


void fill_ref_noind(RefGenome& ref,
                    std::string fasta_file,
                    const bool& cut_names,
                    const bool& remove_soft_mask);


void fill_ref_ind(RefGenome& ref,
                  std::string fasta_file,
                  std::string fai_file,
                  const bool& remove_soft_mask);

#endif
