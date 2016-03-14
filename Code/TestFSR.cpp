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

  return 1;
}
