#include "FSR.h"
#include "TMath.h"
#include <iomanip>      // std::setw
#include <algorithm>    // std::min
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <sstream>
#include <fstream>
//-------------------------------------------------
// Particle
//-------------------------------------------------
Particle::Particle(ParticleType type, double t, double x) : 
  m_type(type)
  , m_t(t)
  , m_x(x)
{}


//-------------------------------------------------
// FSR
//-------------------------------------------------
FSR::FSR(double t0): 
  m_rand(new TRandom3(3.1415))
  , m_t0(t0)
  , m_alpha_s(1./137.) //FIXME: Value
  , m_precision(1e-5)
  , h_t_q(new TH1F("h_t_q","quarks;t;", 100,0,1))
  , h_x_q(new TH1F("h_x_q","quarks;x;", 100,0,1))
  , h_t_g(new TH1F("h_t_g","gluons;t;", 100,0,1))
  , h_x_g(new TH1F("h_x_g","gluons;x;", 100,0,1))
  , h_rnd(new TH1F("h_rnd","random numbers;r",100,0,1))
  , h_radMode(new TH1F("h_radMode",";g->gg, g->qq, q->qg",3, -0.5,2.5))
{  
  DEBUG_MSG("FSR(" << t0 << ")");
}
//-------------------------------------------------
FSR::~FSR()
{
  delete m_rand;
}
//-------------------------------------------------
void FSR::MakeJet(Particle p_in, vector< Particle > &jet)
{
  DEBUG_MSG("FSR::MakeJet");
  m_debugchain.push_back(p_in);    
  if(p_in.GetType() == gluon){
    h_t_g->Fill(p_in.GetT());
    h_x_g->Fill(p_in.GetX());
  }
  else{
    h_t_q->Fill(p_in.GetT());
    h_x_q->Fill(p_in.GetX());
  }

  Particle p_out1, p_out2;
  if(Radiate(p_in, p_out1, p_out2)){
    DEBUG_MSG("; Particle (type, t, x) = (" 
        << p_in.GetType() << ", " << p_in.GetT() << ", " << p_in.GetX() << ") --> radiated --> ("
        << p_out1.GetType() << ", " << p_out1.GetT() << ", " << p_out1.GetX() << ") and ("
        << p_out2.GetType() << ", " << p_out2.GetT() << ", " << p_out2.GetX() << ")");

    MakeJet(p_out1, jet); // call MakeJet recoursively for daughter particles
    //MakeJet(p_out2, jet);
  }
  else
    jet.push_back(p_in);
  return;
}
//-------------------------------------------------
bool FSR::Radiate(Particle p_in, Particle &p_out1, Particle &p_out2)
{
  DEBUG_MSG("FSR::Radiate");

  double t_out1(0), x_out1(0), t_out2(0), x_out2(0);
  ParticleType type_out1(undefined), type_out2(undefined);

  double t_in = p_in.GetT();
  double x_in = p_in.GetX();

  double r_x = m_rand->Uniform(0,1); //TODO: what happens if r_x == 0? What does it mean physically?
  h_rnd->Fill(r_x);

  bool radiated = false;
  if( p_in.GetType() == gluon ){ // is a gluon
    DEBUG_MSG("gluon");
    double r_t_gg = m_rand->Uniform(0,1);
    double r_t_qg = m_rand->Uniform(0,1);
    h_rnd->Fill(r_t_gg);
    h_rnd->Fill(r_t_qg);
 
    double t_out1_gg = GetTFromDelta_gg(t_in, r_t_gg);
    double t_out1_qg = GetTFromDelta_qg(t_in, r_t_qg);
 
//    if(t_out1_gg < t_out1_qg){ // radiate gluons before quarks
      DEBUG_MSG("\t--> radiating 2 gluons");
      h_radMode->Fill(0);
      t_out1 = t_out1_gg;
      x_out1 = GetXFromP_gg(x_in, IntP_gg(0,1) * r_x);
      type_out1 = gluon;
      type_out2 = gluon;
/*    }
    else{ // radiate quarks before gluons
      DEBUG_MSG("\t--> radiating 2 quarks");
      h_radMode->Fill(1);
      t_out1 = t_out1_qg;
      x_out1 = GetXFromP_qg(x_in, IntP_qg(0,1) * r_x);
      type_out1 = quark;
      type_out2 = quark;
    }*/
  }
  else if (p_in.GetType() == quark) { // is a quark
    DEBUG_MSG("quark");
    DEBUG_MSG("\t--> radiating quark and gluon");
    h_radMode->Fill(2);
    
    double r_t_qq = m_rand->Uniform(0,1);
    
    h_rnd->Fill(r_t_qq);
    DEBUG_MSG("r_t_qq " << r_t_qq << ", r_x " << r_x);
    
    t_out1 = GetTFromDelta_qq(t_in, r_t_qq);
    
    
    x_out1 = GetXFromP_qq(x_in, IntP_qq(0,1) * r_x);
    type_out1 = quark;
    type_out2 = gluon;

    DEBUG_MSG("t_out1 " << t_out1 << " x_out1 " << x_out1 );
  }

  t_out2 = TMath::Sqrt(t_in*t_in - t_out1 *t_out1);
  x_out2 = x_in - x_out1;

  if(t_out1 > m_t0 && t_out2 > m_t0 && x_out1 > m_precision && x_out2 > m_precision)
    radiated = true;

  if(radiated){
    p_out1.SetT(t_out1);
    p_out1.SetX(x_out1);
    p_out1.SetType(type_out1);
    p_out2.SetT(t_out2);
    p_out2.SetX(x_out2); // FIXME: check this
    p_out2.SetType(type_out2);

    if(t_out1 > t_in || t_out1 < 0){
      cout << "WARNING: t_out1 > t_in || t_out < 0; t_out1 = " << t_out1 << " t_in = " << t_in << endl; 
      radiated = false;
    }
    if(x_out1 > x_in || x_out1 < 0){
      cout << "WARNING: x_out1 > x_in || x_out < 0; x_out1 = " << x_out1 << " x_in = " << t_in << endl; 
      radiated = false;
    }
  }
  return radiated;
}
//-------------------------------------------------
// Equation (2.166), page 157, g->gg part
double FSR::Delta_gg(double t0, double t1){
  DEBUG_MSG("FSR::Delta_gg");
  double intP_gg = IntP_gg(0,1);
  double retVal = TMath::Exp( -intP_gg * TMath::Log(t1/t0));

  DEBUG_MSG("\t--> " << retVal);
  return retVal;
}
//-------------------------------------------------
// Equation (2.166), page 157, g->qq part
double FSR::Delta_qg(double t0, double t1){
  DEBUG_MSG("FSR::Delta_qg");
  double intP_qg = IntP_qg(0,1);
  double retVal =  TMath::Exp(-intP_qg * TMath::Log(t1/t0));

  DEBUG_MSG("\t--> "<< retVal);
  return retVal;
}
//-------------------------------------------------
// Equation (2.166), page 157, quark part
double FSR::Delta_qq(double t0, double t1)
{
  DEBUG_MSG("FSR::Delta_qq");
  double intP_qq = IntP_qq(0,1);
  double retVal = TMath::Exp( -intP_qq * TMath::Log(t1/t0) );

  DEBUG_MSG("\t--> " << retVal);
  return retVal;
}
//-------------------------------------------------
// Equation (2.105), page 140
double FSR::P_gg(double z)
{
  double C_A = 3; // FIXME: value
  // numerical cutoff:
  if(z<m_integrationCutoff)
    z = m_integrationCutoff;
  else if(z>1-m_integrationCutoff)
    z = 1-m_integrationCutoff;

  double retVal = C_A * ( z/(1-z) + (1-z)/z + z*(1-z) );
  
  return retVal;
}
//-------------------------------------------------
// Equation (2.114), page 142
double FSR::P_qg(double z)
{
  double T_R = 1./2.; // FIXME: value
  double retVal = T_R * ( z*z + (1-z)*(1-z) );

  return retVal;
}
//-------------------------------------------------
// Equation (2.118), page 143
double FSR::P_qq(double z)
{
  double C_F = 4./3.; // FIXME: value
  // numerical cutoff:
  if(z > (1.-m_integrationCutoff))
    z = 1.-m_integrationCutoff;
  double retVal = C_F * (1+z*z) / (1-z);

  return retVal;
}
//-------------------------------------------------
double FSR::Integrate(double (*func)(double), double z0, double z1)
{
  int n = (z1 - z0) / m_precision;
  double z(0), integral(0);

  z1 = z0+m_precision;
  for(int i=0; i<n; i++){
    integral += 0.5 * (m_precision) * ( (*func)(z0+m_precision) + (*func)(z0) );
    z0+=m_precision;
  }
  return /*m_precision * */ integral;
}

double FSR::IntP_gg(double z0, double z1)
{
  DEBUG_MSG("FSR::IntP_gg");
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_gg, z0, z1);
  DEBUG_MSG("\t--> " << retVal);
  return retVal;
}
//-------------------------------------------------
double FSR::IntP_qg(double z0, double z1)
{
  DEBUG_MSG("FSR::IntP_qg");
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_qg, z0, z1);
  DEBUG_MSG("\t--> " << retVal);
  return retVal;
}
//-------------------------------------------------
double FSR::IntP_qq(double z0, double z1)
{
  DEBUG_MSG("FSR::IntP_qq");
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_qq, z0, z1);
  DEBUG_MSG("\t--> " << retVal);

  return retVal;
}
//-------------------------------------------------
bool FSR::CanRadiate(Particle p)
{
  DEBUG_MSG( "FSR::CanRadiate" );
  if(p.GetT() > m_t0 && p.GetX() > m_precision) // FIXME: in principle x > 0, but how precise can we compute it?
    return true;
  else 
    return false;
}
//-------------------------------------------------
// Equation (2.174), page 159
double FSR::GetTFromDelta_gg(double t_in, double c)
{
  DEBUG_MSG("FSR::GetTFromDelta_gg");
  double intP_gg = IntP_gg(0,1);
  double t_out(t_in);
  double arg = intP_gg * TMath::Log(t_in/t_out);

  while(TMath::Exp(-arg) > c){
    arg = intP_gg * TMath::Log(t_in/t_out);
    t_out -= m_precision;
  }
  DEBUG_MSG("\t--> " << t_out );
  return t_out;
}
//-------------------------------------------------
// Equation (2.174), page 159
double FSR::GetTFromDelta_qg(double t_in, double c)
{
  DEBUG_MSG("FSR::GetTFromDelta_qg");
  double intP_qg = IntP_qg(0,1);
  double t_out(t_in);
  double arg = intP_qg * TMath::Log(t_in/t_out);


  while(TMath::Exp(-arg) > c && t_out >= m_t0){
    arg = intP_qg * TMath::Log(t_in/t_out);
    t_out -= m_precision;
  }
  DEBUG_MSG("\t--> " << t_out );
  return t_out;
}
//-------------------------------------------------
// Equation (2.174), page 159
double FSR::GetTFromDelta_qq(double t_in, double c){
  DEBUG_MSG("FSR::GetTFromDelta_qq");
  double intP_qq = IntP_qq(0,1);
  double t_out(t_in);
  double arg = intP_qq * TMath::Log(t_in/t_out);

  intP_qq = 10;

  while(TMath::Exp(-arg) > c && t_out >= m_t0){
   
    
    arg = intP_qq * TMath::Log(t_in/t_out);
    
    
    t_out -= m_precision;
  }
  DEBUG_MSG("\t--> " << t_out );
  return t_out;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_gg(double x1, double c)
{
  DEBUG_MSG( "FSR::GetXFromP_gg" );
  double integral(0);
  double z0(0);
  double z1 = z0+m_precision;
  while(integral < c){
    integral += Integrate(&P_gg, z0, z0+m_precision);
    z0 += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  double x_out = (2*z0+m_precision)/2. * x1;
  DEBUG_MSG("\t-->"<<x_out);
  return x_out;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_qg(double x1, double c) 
{
  DEBUG_MSG( "FSR::GetXFromP_qg" );
  double integral(0);
  double z0(0);
  while(integral < c){ 
    integral += Integrate(&P_qg, z0, z0+m_precision);
    z0 += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  double x_out = (2*z0+m_precision)/2. * x1;
  DEBUG_MSG("\t-->"<<x_out);
  return x_out;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_qq(double x1, double c) 
{
  DEBUG_MSG( "FSR::GetXFromP_qq" );
  double integral(0);
  double z0(0);
    while(integral < c){ 
      integral += Integrate(&P_qq, z0, z0 + m_precision);
      z0 += m_precision;
    }
  // upper integration bound is something x0 == x2/x1
  double x_out = (2*z0+m_precision)/2. * x1;
  DEBUG_MSG("\t-->"<<x_out);
  return x_out;
}
//-------------------------------------------------
void FSR::DrawTXPlot(char* pdf){
  TCanvas c("c");
  c.Divide(2,2);
  vector< double > t;
  vector< double > x;
  int n = min((int)m_debugchain.size(),100);
  for(int i=0; i<n; i++){
    t.push_back(m_debugchain[i].GetT());
    x.push_back(m_debugchain[i].GetX());
  }
  TGraph g(n, &t[0], &x[0]);
  if(m_debugchain[0].GetType() == quark)
    g.SetTitle("Initial quark");
  if(m_debugchain[0].GetType() == gluon)
    g.SetTitle("Initial gluon");
  c.cd(1);
  g.Draw("A*");
  m_debugchain.clear();
  
  c.cd(2);
  h_radMode->Draw();

  c.cd(3);
  h_t_q->Draw();
  h_t_g->SetLineColor(2);
  h_t_g->Draw("same");
  gPad->BuildLegend()->Draw();
  c.cd(4);
  h_x_q->Draw();
  h_x_g->SetLineColor(2);
  h_x_g->Draw("same");
  c.Print(pdf);

}

// Load and save events //

void FSR::save_events(std::string filename, const std::vector<event>& events)
{
  TFile* output = new TFile(filename.c_str(), "RECREATE");
  output->cd();
  TTree* tree = new TTree("Events", "Events");
  event meas;  
  tree->Branch("weight",&meas.weight,"weight/D");
  tree->Branch("part","std::vector<Particle>",&meas.part);
  for (unsigned int i=0; i< events.size(); i++)
    {
      meas = events.at(i);
      tree->Fill();
    }
  tree->Write();
  output->Write();
  output->Close();
  delete output;
}


std::vector<event> FSR::load_events(std::string filename, int neventsmax)
{
  ifstream input;
  input.open(filename.c_str());
  std::cout << "Reading input file: " << filename << std::endl;
  event meas;  
  std::vector<event> result;
  
  string line;
  int n_events = 0;
  //Read in events from text file
  if (input.is_open() ){
      //loop over all lines
      while ( getline (input,line) ){

          //check number of read events (=lines) and break if requested maximum reached
          n_events++;
          if (n_events > neventsmax) break;

          DEBUG_MSG("Line read in:" << line);

          //set number of gluons/quarks (0/2 per default)
          meas.nGluons = 0;
          meas.nQuarks = 2;

          std::vector<double> v;
          std::string s;
          istringstream linestream;
          linestream.str( line );
          //decompose line string into required properties and push into double vector
          while (getline (linestream, s, ' ')){
              v.push_back((double)atof(s.c_str()));
          }

          //finally, take vector and fill event properties (super hardcoded, thanks Johann et al.)
          //event weight
          DEBUG_MSG("Event" << n_events << ": weight=" << v.at(0) );
          meas.weight = v.at(0);
          //first quark
          double x1 = TMath::Sqrt( v.at(2)*v.at(2) + v.at(3)*v.at(3) + v.at(4)*v.at(4) );
          double t1 = TMath::Sqrt( v.at(2)*v.at(2) + v.at(3)*v.at(3) ) / x1; //FIXME: using PT/P, correct?
          double E1 = v.at(1);
          Particle quark1 = Particle(quark, t1, x1);
          meas.part.push_back(quark1);
          //second quark
          double x2 = TMath::Sqrt( v.at(6)*v.at(6) + v.at(7)*v.at(7) + v.at(8)*v.at(8) );
          double t2 = TMath::Sqrt( v.at(6)*v.at(6) + v.at(7)*v.at(7) ) / x2; //FIXME: using PT/P, correct?
          double E2 = v.at(5);
          Particle quark2 = Particle(quark, t2, x2);
          meas.part.push_back(quark2);

          result.push_back(meas);

      }

  }
  input.close();

  return result;
}
