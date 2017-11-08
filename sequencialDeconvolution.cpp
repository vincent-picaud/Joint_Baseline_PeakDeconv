// The "more usual approach":
// - smoothing
// - SNIP (baseline removal)
// - l1-deconvolution
//
#include <csv.h>
#include <cxxopts.hpp>

#include <fstream>
#include <sstream>

namespace JointDeconv
{
};  // TODO: to remove

using namespace JointDeconv;

int main(int argc, char* argv[])
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
  cxxopts::Options options(argv[0],
                           "A more classical approach "
                           "smoothing->SNIP->deconvolution (for comparison "
                           "purpose only), contact vincent.picaud@cea.fr");

  // Default parameter values
  //
  const int half_width_SG = 19;
  double SG[2 * half_width_SG + 1];  // TODO

  options.add_options()
      //
      ("i,input",
       "Input file (two columns X,Y)",
       cxxopts::value<std::string>(),
       "FILE")
      //
      ("help", "Print help");

  options.parse_positional("input");

  // CAVEAT: after options.parse(argc, argv); argc is set=1
  // -> we have to check its value before
  //
  bool show_help = argc == 1;

  options.parse(argc, argv);

  show_help |= options.count("help") > 0;

  if (show_help)
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  // Sanity check
  //
  if (options.count("input") == 0)
  {
    std::cout << "#Error: missing input FILE" << std::endl;
    exit(-1);
  }
  // ...

  // Generates output filename
  //
  const std::string output_filename((options.count("output") == 0)
                                        ? options["input"].as<std::string>() +
                                              ".out"
                                        : options["output"].as<std::string>());
}
