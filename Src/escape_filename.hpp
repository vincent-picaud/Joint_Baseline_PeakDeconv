#pragma once

#include <string>

namespace JointDeconv
{
  // For filename as gnuplot title, need to escape _
  std::string escape_filename(const std::string source)
  {
    std::string escaped;
    for (auto c : source)
    {
      if (c == '_')
      {
        escaped += "\\";
      }
      escaped += c;
    }
    return escaped;
  }
}  // namespace JointDeconv
