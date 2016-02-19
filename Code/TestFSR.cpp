#include "FSR.h"

#include "iostream"
#include "math.h"
#include <iomanip>      // std::setw
#include <TH1.h>

using namespace std;

int main()
{
 

  // FSR Test
  FSR test(0.001);
  Particle p_quark(quark, 1, 1);
  Particle p_gluon(gluon, 1, 1);
  vector< Particle > jet_out;
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf(");
  jet_out.clear();
  test.MakeJet(p_gluon, jet_out);
  test.DrawTXPlot("TestFSR.pdf)");


  return 1;
}
