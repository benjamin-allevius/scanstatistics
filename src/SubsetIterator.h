#ifndef SUBSETITERATOR_H
#define SUBSETITERATOR_H

#include <iterator>
#include <memory>
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

class SubsetIterator : public std::iterator<std::input_iterator_tag, arma::uvec>
{
public:
  SubsetIterator(std::shared_ptr<arma::uvec> subset_ptr);
  
private: 
  
  
};


#endif
