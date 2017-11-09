#pragma once

#include "vector.hpp"

namespace JointDeconv
{
  void sgFilter(const Vector& spectrum, Vector& filtered);
}  // JointDeconv
