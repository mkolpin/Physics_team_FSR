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

 vector<event> events = test.load_events("events.dat",-1);
  vector<event> events_out;
  for (auto evt : events){
      vector<Particle> out_parts;
      //event event1 = events.at(1);
      for (auto part : evt.part){
          //Particle p_quark = evt.part.at(1);
          Particle p_quark = part;
          Particle p_gluon(gluon, 1, 1);
          vector< Particle > jet_out;
          test.MakeJet(p_quark, jet_out);
          test.DrawTXPlot("TestFSR.pdf(");
          jet_out.clear();
          test.MakeJet(p_quark, jet_out);
          test.DrawTXPlot("TestFSR.pdf)");
          out_parts.insert(out_parts.end(), jet_out.begin(), jet_out.end());
      }
      event evt_out;
      evt_out.weight = evt.weight;
      evt_out.part = out_parts;
      events_out.push_back(evt_out);
  }

  //save to output root file
  test.save_events("test.root",events_out);





  c.Print("TestFSR.pdf]");
  return 1;
}
