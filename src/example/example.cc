#include <fstream>
#include <iostream>

#include "cgb.hh"

int main(int argc, char **argv) {
  if (argc < 3 || argc > 6) {
    std::cerr << "Usage: " << argv[0]
              << " input.cgb output.obj [resolution] [levels]"
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

  {
    std::ifstream f(argv[1]);
    if (!f.is_open() || !cgb.loadOptions(f) || !cgb.loadControlPoints(f)) {
      std::cerr << "Cannot open file: " << argv[1] << std::endl;
      return 2;
    }
  }

#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif

  cgb.evaluate(resolution).writeOBJ(argv[2]);

  return 0;
}
