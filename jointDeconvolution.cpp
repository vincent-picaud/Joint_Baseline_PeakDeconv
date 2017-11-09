#include "jointDeconv.hpp"

#include <csv.h>
#include <cxxopts.hpp>

#include <fstream>
#include <sstream>

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
                           "A joint baseline removal and deconvolution "
                           "algorithm, contact vincent.picaud@cea.fr");

  // Default parameter values
  //
  double sigma_left = 10;
  double sigma_right = 10;
  double peakMinHeight = 0.01;
  double yb_left;
  double yb_right;

  double lambda_1 = 0.1;
  double lambda_2 = 0.00001;
  double mu = 500;
  double eps = 1e-4;
  Size_t max_iter = 5000;
  bool gnuplot = false;

  options.add_options()
      //
      ("i,input",
       "Input file (two columns X,Y)",
       cxxopts::value<std::string>(),
       "FILE")
      //
      ("o,output",
       "Output file",
       cxxopts::value<std::string>()->default_value("$(FILE).out"),
       "OUTPUT FILE")
      //
      ("sigma_left",
       "Peak shape factor (>0)",
       cxxopts::value<double>(sigma_left),
       to_string(sigma_left))
      //
      ("sigma_right",
       "Peak shape factor (>0)",
       cxxopts::value<double>(sigma_right),
       to_string(sigma_right))
      //
      ("yb_left",
       "Left baseline value (if not defined use y[0])",
       cxxopts::value<double>(yb_left),
       "y[0]")
      //
      ("yb_right",
       "Right baseline value (if not defined use y[n-1])",
       cxxopts::value<double>(yb_right),
       "y[n-1]")
      //
      ("peakMinHeight",
       "Minimal height to accept peak (>=0)",
       cxxopts::value<double>(peakMinHeight),
       to_string(peakMinHeight))
      //
      ("lambda_1",
       "lambda_1 penalty term (>=0)",
       cxxopts::value<double>(lambda_1),
       to_string(lambda_1))
      //
      ("lambda_2",
       "lambda_2 penalty term (>=0)",
       cxxopts::value<double>(lambda_2),
       to_string(lambda_2))
      //
      ("mu", "mu penalty term (>0)", cxxopts::value<double>(mu), to_string(mu))
      //
      ("eps", "eps goal (>=0)", cxxopts::value<double>(eps), to_string(eps))
      //
      ("max_iter",
       "maximum number of iterations (>0)",
       cxxopts::value<Size_t>(max_iter),
       to_string(max_iter))
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
    return EXIT_FAILURE;
  }
  if (sigma_left <= 0)
  {
    std::cout << "#Error: sigma_left= " << sigma_left
              << " is not a positive number" << std::endl;
    return EXIT_FAILURE;
  }
  if (sigma_right <= 0)
  {
    std::cout << "#Error: sigma_right= " << sigma_right
              << " is not a positive number" << std::endl;
    return EXIT_FAILURE;
  }
  if (peakMinHeight < 0)
  {
    std::cout << "#Error: peakMinHeight= " << peakMinHeight
              << " is not a nonegative number" << std::endl;
    return EXIT_FAILURE;
  }
  if (lambda_1 < 0)
  {
    std::cout << "#Error: lambda_1= " << lambda_1
              << " is not a nonnegative number" << std::endl;
    return EXIT_FAILURE;
  }
  if (lambda_2 < 0)
  {
    std::cout << "#Error: lambda_2= " << lambda_2
              << " is not a nonnegative number" << std::endl;
    return EXIT_FAILURE;
  }
  if (mu <= 0)
  {
    std::cout << "#Error: mu= " << mu << " is not a positive number"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (eps < 0)
  {
    std::cout << "#Error: eps= " << eps << " is not a nonnegative number"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (max_iter <= 0)
  {
    std::cout << "#Error: max_iter= " << max_iter << " is not a positive number"
              << std::endl;
    return EXIT_FAILURE;
  }

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
    return EXIT_FAILURE;
  }

  ///////////////////////////////////////
  // Call the deconvolution subroutine //
  ///////////////////////////////////////

  const Size_t n = data.size();

  if (n <= 2)
  {
    std::cerr << "#Error: spectrum must have at least 2 points" << std::endl;
    return EXIT_FAILURE;
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

  JointDeconv_InputParameters inputParameters;
  inputParameters.lambda_1 = lambda_1;
  inputParameters.lambda_2 = lambda_2;
  inputParameters.mu = mu;
  inputParameters.solver_inputParameters.iter_max = max_iter;
  inputParameters.solver_inputParameters.eps_goal = eps;
  jointDeconv_GaussianPeaks(
      x,
      y,
      (options.count("yb_left") > 0) ? yb_left : y[0],
      (options.count("yb_right") > 0) ? yb_right : y[n - 1],
      peakMinHeight,
      sigma_left,
      sigma_right,
      deconvolvedPeak,
      convolvedPeak,
      baseline,
      inputParameters);

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

  return EXIT_SUCCESS;
}
