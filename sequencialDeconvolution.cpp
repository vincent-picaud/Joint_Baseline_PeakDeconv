// The "more usual approach":
// - smoothing
// - SNIP (baseline removal)
// - l1-deconvolution
//
#include <csv.h>
#include <cxxopts.hpp>

#include <fstream>
#include <sstream>

#include "snip.hpp"

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

  bool gnuplot = false;

  options.add_options()
      //
      ("i,input",
       "Input file (two columns X,Y)",
       cxxopts::value<std::string>(),
       "FILE")
      //
      ("p,gnuplot", "Gnuplot script", cxxopts::value<bool>(gnuplot))
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

  ///////////////////
  // Read CSV file //
  ///////////////////

  std::vector<std::pair<double, double>> data;
  try
  {
    io::CSVReader<2,
                  io::trim_chars<' ', '\t'>,
                  io::no_quote_escape<','>,
                  io::throw_on_overflow,
                  io::single_line_comment<'#'>>
        in(options["input"].as<std::string>());
    double X_i, Y_i;

    while (in.read_row(X_i, Y_i))
    {
      data.push_back(std::make_pair(X_i, Y_i));
    }
  }
  catch (std::exception& e)
  {
    std::cerr << "#Error: " << e.what() << std::endl;
    exit(-1);
  }

  ///////////////////////////////////////
  // Call the deconvolution subroutine //
  ///////////////////////////////////////

  const Size_t n = data.size();

  if (n <= 2)
  {
    std::cerr << "#Error: spectrum must have at least 2 points" << std::endl;
    exit(-1);
  }

  Vector x(n);
  Vector y(n);
  Vector deconvolvedPeak(n);
  Vector convolvedPeak(n);
  Vector baseline(n);

  for (Index_t i = 0; i < n; ++i)
  {
    x[i] = data[i].first;
    y[i] = data[i].second;
  }

  data.clear();

  //////////////////
  // Computation TODO //
  //////////////////

  // SNIP
  snip(y, baseline, 20);

  //////////////////
  // Write output //
  //////////////////

  std::ofstream output(output_filename);
  output << "# Generated by: " << command_line << std::endl;
  for (Index_t i = 0; i < n; ++i)
  {
    output << x[i] << ", " << y[i] << ", " << deconvolvedPeak[i] << ", "
           << baseline[i] << ", " << convolvedPeak[i] << ", "
           << baseline[i] + convolvedPeak[i] << std::endl;
  }

  output.close();

  // Gnuplot
  // load "test.csv.out.gnuplot"
  if (gnuplot)
  {
    std::ofstream output_gnuplot(output_filename + ".gnuplot");

    output_gnuplot << "set datafile separator ','" << std::endl;
    output_gnuplot << "plot \"" << output_filename
                   << "\" u 1 : 2 w l t \"Y raw\" lc \"blue\"" << std::endl;
    output_gnuplot << "replot \"" << output_filename
                   << "\" u 1 : 6 w l t \"reconstructed\" lw 2 lc \"grey30\""
                   << std::endl;
    output_gnuplot << "replot \"" << output_filename
                   << "\" u 1 : 4 w l t \"baseline\" lw 2 lc \"grey70\""
                   << std::endl;
    output_gnuplot << "replot \"" << output_filename
                   << "\" u 1 : 5 w l t \"convolved peaks\" lc \"green\""
                   << std::endl;
    output_gnuplot << "replot \"" << output_filename
                   << "\" u 1 : 3 w i t \"peaks\" lw 2 lc \"red\"" << std::endl;
    output_gnuplot << "set terminal eps" << std::endl;
    output_gnuplot << "set output \"" << output_filename + ".eps"
                   << "\"" << std::endl;
    output_gnuplot << "replot" << std::endl;
    output_gnuplot << "set terminal png" << std::endl;
    output_gnuplot << "set output \"" << output_filename + ".png"
                   << "\"" << std::endl;
    output_gnuplot << "replot" << std::endl;
    output_gnuplot << "set terminal qt" << std::endl;
    output_gnuplot << "set output" << std::endl;
  }
}