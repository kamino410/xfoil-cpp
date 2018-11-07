#include <iostream>
#include <sstream>

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

int main() {
  auto naca = new XFoil();
  naca->naca4(4412, 100 / 2);

  double xb[604], yb[604], x[604], y[604], nx[604], ny[604];
  for (int i = 0; i < naca->nb; i++) {
    xb[i] = naca->xb[i + 1];
    yb[i] = naca->yb[i + 1];
    x[i] = naca->xb[i + 1];
    y[i] = naca->yb[i + 1];
  }

  XFoil *foil = new XFoil();

  std::stringstream ss;
  if (!foil->initXFoilGeometry(naca->nb, x, y, nx, ny)) {
    std::cout << "Initialization error!" << std::endl;
    return 0;
  }
  if (!foil->initXFoilAnalysis(100000, 0, 0.0, 9.0, 1.0, 1.0, 1, 1, true, ss)) {
    std::cout << "Initialization error!" << std::endl;
    return 0;
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
      getchar();
      return 0;
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
