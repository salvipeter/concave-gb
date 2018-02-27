#include <fstream>
#include <iostream>

#include "cgb.hh"

int main(int argc, char **argv) {
  if (argc < 3 || argc > 7) {
    std::cerr << "Usage: " << argv[0]
              << " input.cgb output.obj [resolution] [levels] [biharmonic] [fill-concave]"
              << std::endl;
    return 1;
  }

  CGB::ConcaveGB cgb;

  double resolution = 0.001;
  
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
      cgb.setBiharmonic(true);
  }

  if (argc > 6) {
    if (std::string(argv[6]) == "true")
      cgb.setFillConcaveCorners(true);
  }

#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif
  std::cout << "Input: " << argv[1] << std::endl;
  std::cout << "Output: " << argv[2] << std::endl;
  std::cout << "Resolution: " << resolution << std::endl;
  std::cout << "Param. level: " << (argc > 4 ? argv[4] : "9") << std::endl;
  std::cout << "Biharmonic: " << (argc > 5 && std::string(argv[5]) == "true" ? "true" : "false")
            << std::endl;
  std::cout << "Fill corners: " << (argc > 6 && std::string(argv[6]) == "true" ? "true" : "false")
            << std::endl;

  {
    std::ifstream f(argv[1]);
    if (!f.is_open() || !cgb.loadOptions(f) || !cgb.loadControlPoints(f)) {
      std::cerr << "Cannot open file: " << argv[1] << std::endl;
      return 2;
    }
  }

  cgb.evaluate(resolution).writeOBJ(argv[2]);

  return 0;
}
