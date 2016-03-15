#include "FSR.h"

#include "iostream"
#include "math.h"
#include <iomanip>      // std::setw
#include <TH1.h>
#include <TCanvas.h>
#include "TStyle.h"

using namespace std;
  
int main()
{
  gStyle->SetOptStat(111111);

  TCanvas c("c");
  c.Print("TestFSR.pdf[");
 

  // FSR Test
  double t_in = 100;
  double t0 = 1;
  FSR test(1);
  Particle p_quark(quark, t_in, 1);
  Particle p_gluon(gluon, t_in, 1);
  vector< Particle > jet_out;

  test.DebugPlots("TestFSR.pdf", t_in);
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf");
  jet_out.clear();
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf");
  jet_out.clear();
  cout << "\n=================== NEW JET =======================\n";
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf");

 
  TH1F* h_jet_t = new TH1F("h_jet_t", ";t",100,0,t_in+0.1);
  TH1F* h_jet_x = new TH1F("h_jet_x", ";x",100,0,1.1);
  TH1F* h_jet_mult = new TH1F("h_jet_mult", ";#jets",10,-.5,9.5);
  for(int i=0; i<100; i++){
    if(i%10 == 0) cout << "=== " << i << " ===\n";
    jet_out.clear();
    test.MakeJet(p_quark, jet_out);

    for(auto p: jet_out){
      h_jet_t->Fill(p.GetT());
      h_jet_x->Fill(p.GetX());
    }
    h_jet_mult->Fill(jet_out.size());
  }
  c.Divide(2,2);
  c.cd(1);
  h_jet_t->Draw();
  c.cd(2);
  h_jet_x->Draw();
  c.cd(3);
  h_jet_mult->Draw();
  c.Print("TestFSR.pdf");
  
  test.DrawTXPlot("TestFSR.pdf");


  c.Print("TestFSR.pdf]");
  return 1;
}
