#include "syntheticExample.hpp"

#include <csv.h>
#include <cxxopts.hpp>

#include <fstream>
#include <sstream>

using namespace JointDeconv;

int
main(int argc, char* argv[])
{
  //////////////////////
  // Argument parsing //
  //////////////////////

  const auto to_string = [](const auto& any) {
    std::stringstream converter;
    converter << any;
    return std::string(converter.str());
  };

  std::string command_line;
  for (int i = 0; i < argc; i++)
  {
    command_line += argv[i];
    command_line += ' ';
  }

  cxxopts::Options options(argv[0], "Generates some synthetic data");

  // Default parameter values
  //
  int    c        = +1;
  double noiseSTD = 0.2;

  options.add_options()
      //
      ("c,concavity", "Baseline concavity, an integer in {-1,0,+1}", cxxopts::value<int>(c), to_string(c))
      //
      ("n,noise", "Noise standard deviation (>0)", cxxopts::value<double>(noiseSTD), to_string(noiseSTD))
      //
      ("help", "Print help");

  // CAVEAT: after options.parse(argc, argv); argc is set=1
  // -> we have to check its value before
  //
  options.parse(argc, argv);

  if (options.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  if ((c != -1) && (c != 0) && (c != +1))
  {
    std::cerr << "\nError: c value " << c << " is not in {-1,0,+1}" << std::endl;
    exit(-1);
  }

  if (noiseSTD < 0)
  {
    std::cout << "#Error: noiseSTD= " << noiseSTD << " is not a nonnegative number" << std::endl;
    exit(-1);
  }

  const auto data = create_syntheticExample(c, noiseSTD);
  const auto n    = data.y.size();

  std::cout << "# " << command_line << std::endl;

  // To avoid complication use the true baseline value at boundaries
  std::cout << 0 << "," << data.y_baseline[0] << std::endl;
  for (Index_t i = 1; i + 1 < n; i++)
  {
    std::cout << i << "," << data.y[i] << std::endl;
  }
  std::cout << n - 1 << "," << data.y_baseline[n - 1] << std::endl;
}
