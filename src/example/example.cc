#include <fstream>
#include <iostream>

#include "cgb.hh"

#ifdef WIN32
extern "C" { FILE __iob_func[3] = { *stdin, *stdout, *stderr }; } // GSL is compiled with an old VSs
#endif

int main(int argc, char **argv) {
  if (argc < 3 || argc > 8) {
    std::cerr << "Usage: " << argv[0]
              << " input.cgb output.obj"
              << " [resolution] [levels] [maxent] [fill-concave] [concave-weight]"
              << std::endl;
    return 1;
  }

  CGB::ConcaveGB cgb;

  double resolution = 0.001;
  double concave_weight = 1.0;
  
  if (argc > 3)
    resolution = std::stod(argv[3]);

  if (argc > 4) {
    unsigned int levels = std::stoi(argv[4]);
    if (levels <= 0) {
      std::cerr << "Parameterization level should be positive." << std::endl;
      return 1;
    }
    cgb.setParamLevels(levels);
  }

  if (argc > 5) {
    if (std::string(argv[5]) == "true")
      cgb.setMaxEnt(true);
  }

  if (argc > 6) {
    if (std::string(argv[6]) == "true")
      cgb.setFillConcaveCorners(true);
  }

  if (argc > 7) {
    concave_weight = std::stod(argv[7]);
  }
  cgb.setConcaveWeight(concave_weight);

#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif
  std::cout << "Input: " << argv[1] << std::endl;
  std::cout << "Output: " << argv[2] << std::endl;
  std::cout << "Resolution: " << resolution << std::endl;
  std::cout << "Param. level: " << (argc > 4 ? argv[4] : "9") << std::endl;
  std::cout << "Max. entropy: " << (argc > 5 && std::string(argv[5]) == "true" ? "true" : "false")
            << std::endl;
  std::cout << "Fill corners: " << (argc > 6 && std::string(argv[6]) == "true" ? "true" : "false")
            << std::endl;
  std::cout << "Concave weight: " << concave_weight << std::endl;

  {
    std::ifstream f(argv[1]);
    if (!f.is_open() || !cgb.loadOptions(f) || !cgb.loadControlPoints(f)) {
      std::cerr << "Cannot open file: " << argv[1] << std::endl;
      return 2;
    }
  }

  if (resolution > 0)
    cgb.evaluate(resolution).writeOBJ(argv[2]);
  else
    cgb.evaluateRegular(-resolution).writeOBJ(argv[2]);

  return 0;
}
