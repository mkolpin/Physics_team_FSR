#include "FSR.h"

#include "iostream"
#include "math.h"
#include <iomanip>      // std::setw
#include <TH1.h>
#include <TCanvas.h>

using namespace std;
  
int main()
{
  TCanvas c("c");
  c.Print("TestFSR.pdf[");
 

  // FSR Test
  FSR test(1e-3);
  Particle p_quark(quark, 1, 1);
  Particle p_gluon(gluon, 1, 1);
  vector< Particle > jet_out;

  test.DebugPlots("TestFSR.pdf");
  
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf");
  jet_out.clear();
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf");
  jet_out.clear();
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_gluon, jet_out);
  test.DrawTXPlot("TestFSR.pdf");

  TH1F* h_jet_t = new TH1F("h_jet_t", "",100,0,1);
  TH1F* h_jet_x = new TH1F("h_jet_x", "",100,0,1);
  for(int i=0; i<1000; i++){
    jet_out.clear();
    test.MakeJet(p_quark, jet_out);

    cout << jet_out.size() << endl;
    for(auto p: jet_out){
      h_jet_t->Fill(p.GetT());
      h_jet_x->Fill(p.GetX());
    }
  }
  c.Divide(2,2);
  c.cd(1);
  h_jet_t->Draw();
  c.cd(2);
  h_jet_x->Draw();
  c.Print("TestFSR.pdf");

  c.Print("TestFSR.pdf]");
  return 1;
}
