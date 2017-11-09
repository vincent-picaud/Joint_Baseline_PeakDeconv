#pragma once

#include "snip.hpp"
#include "vector.hpp"

namespace JointDeconv
{
  void snip(const Vector& spectrum, Vector& baseline, const Size_t halfWindowSize);
}  // JointDeconv
