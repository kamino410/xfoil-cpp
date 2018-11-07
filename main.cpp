#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "XFoil.h"

int m_Iterations = 0;
int s_IterLim = 100;
bool m_bErrors = false;
bool s_bAutoInitBL = true;

bool iterate(XFoil *foil) {
  if (!foil->viscal()) {
    foil->lvconv = false;
    std::cout
        << "CpCalc: local speed too large\nCompressibility corrections invalid"
        << std::endl;
    return false;
  }

  while (m_Iterations < s_IterLim && !foil->lvconv /*&& !s_bCancel*/) {
    if (foil->ViscousIter()) {
      // if (m_x0 && m_y0) {
      //  m_x0->append((double)m_Iterations);
      //  m_y0->append(foil->rmsbl);
      //}
      // if (m_x1 && m_y1) {
      //  m_x1->append((double)m_Iterations);
      //  m_y1->append(foil->rmxbl);
      //}
      m_Iterations++;
    } else
      m_Iterations = s_IterLim;
  }

  if (!foil->ViscalEnd()) {
    foil->lvconv = false;  // point is unconverged

    foil->setBLInitialized(false);
    foil->lipan = false;
    m_bErrors = true;
    return true;  // to exit loop
  }

  if (m_Iterations >= s_IterLim && !foil->lvconv) {
    if (s_bAutoInitBL) {
      foil->setBLInitialized(false);
      foil->lipan = false;
    }
    foil->fcpmin();  // Is it of any use ?
    return true;
  }
  if (!foil->lvconv) {
    m_bErrors = true;
    foil->fcpmin();  // Is it of any use ?
    return false;
  } else {
    // converged at last
    foil->fcpmin();  // Is it of any use ?
    return true;
  }
  return false;
}

int loadDatFile(std::string filename, double x[604], double y[604]) {
  std::ifstream fs(filename);
  if (!fs) {
    std::cout << "Failed to open dat file" << std::endl;
    return -1;
  }

  std::string line;
  std::getline(fs, line);
  std::cout << "Foil name : " << line << std::endl;
  int cnt = 0;
  while (!fs.eof()) {
    std::getline(fs, line);

    line.erase(0, line.find_first_not_of(" \t"));
    int endOfX = line.find_first_of(" \t");
    if (endOfX == -1) continue;

    std::string sx = line.substr(0, endOfX);
    std::string sy = line.substr(endOfX);

    x[cnt] = atof(sx.c_str());
    y[cnt] = atof(sy.c_str());
    cnt++;
  }
  return cnt;
}

int main() {
  double x[604], y[604], nx[604], ny[604];
  int n = 0;

  std::stringstream ss;

  // auto naca = new XFoil();
  // naca->naca4(4412, 100 / 2);
  // for (int i = 0; i < naca->nb; i++) {
  //   x[i] = naca->xb[i + 1];
  //   y[i] = naca->yb[i + 1];
  // }
  // n = naca->nb;

  n = loadDatFile("sample/CLARK_Y.dat", x, y);
  if (n == -1) return 1;

  XFoil *foil = new XFoil();

  if (!foil->initXFoilGeometry(n, x, y, nx, ny)) {
    std::cout << "Initialization error!" << std::endl;
    return 1;
  }
  if (!foil->initXFoilAnalysis(100000, 0, 0.0, 9.0, 1.0, 1.0, 1, 1, true, ss)) {
    std::cout << "Initialization error!" << std::endl;
    return 1;
  }

  for (double alpha = 0; alpha < 15; alpha += 0.5) {
    m_Iterations = 0;

    foil->setBLInitialized(false);
    foil->lipan = false;

    foil->setAlpha(alpha * 3.14159 / 180);
    foil->lalfa = true;
    foil->setQInf(1.0);
    std::cout << "alpha : " << alpha << std::endl;

    if (!foil->specal()) {
      std::cout << "Invalid Analysis Settings" << std::endl;
      return 1;
    }
    foil->lwake = false;
    foil->lvconv = false;

    while (!iterate(foil))
      ;

    // std::cout << ss.str() << std::endl;

    if (foil->lvconv) {
      std::cout << "  converged after " << m_Iterations << " iterations"
                << std::endl;
      std::cout << "  cl : " << foil->cl << ", cd : " << foil->cd
                << ", cm : " << foil->cm << ", xcp : " << foil->xcp
                << std::endl;
    } else {
      std::cout << "  unconverged" << std::endl;
    }
  }

  return 0;
}
