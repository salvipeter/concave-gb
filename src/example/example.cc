#include <iostream>

#include "cgb.hh"

int main(int argc, char **argv) {
  if (argc < 4 || argc > 6) {
    std::cerr << "Usage: " << argv[0]
              << " input.cgb output.obj resolution [central_weight] [levels]"
              << std::endl;
    return 1;
  }

  CGB::ConcaveGB cgb;

  double resolution = std::stod(argv[3]);
  if (resolution <= 0.0) {
    std::cerr << "Resolution should be positive." << std::endl;
    return 1;
  }

  if (argc > 4) {
    std::string type(argv[4]);
    if (type == "original")
      cgb.setCentralWeight(CGB::ConcaveGB::CentralWeight::ORIGINAL);
    else if (type == "zero")
      cgb.setCentralWeight(CGB::ConcaveGB::CentralWeight::ZERO);
    else {
      std::cerr << "Possible values for central_weight: original / zero ." << std::endl;
      return 1;
    }
  }

  if (argc > 5) {
    unsigned int levels = std::stoi(argv[5]);
    if (levels <= 0) {
      std::cerr << "Parameterization level should be positive." << std::endl;
      return 1;
    }
    cgb.setParamLevels(levels);
  }

  if (!cgb.loadControlPoints(argv[1])) {
    std::cerr << "Cannot open file: " << argv[1] << std::endl;
    return 2;
  }

  cgb.evaluate(resolution).writeOBJ(argv[2]);

  return 0;
}
