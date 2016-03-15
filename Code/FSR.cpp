#include "FSR.h"
#include "TMath.h"
#include <iomanip>      // std::setw
#include <algorithm>    // std::min
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
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
 { 
   DEBUG_MSG("FSR(" << t0 << ")");
   h_t_q = new TH1F("h_t_q","quarks;t;", 100,0,plotTMax);
   h_x_q = new TH1F("h_x_q","quarks;x;", 100,0,plotXMax);
   h_t_g = new TH1F("h_t_g","gluons;t;", 100,0,plotTMax);
   h_x_g = new TH1F("h_x_g","gluons;x;", 100,0,plotXMax);
   h_rnd = new TH1F("h_rnd","random numbers;r",100,0,1.);
   h_xIn_xOut = new TH2F("h_xIn_xOut",";x_in;x_out",100,0,plotXMax,100,0,plotXMax);
   h_tIn_tOut = new TH2F("h_tIn_tOut",";t_in;t_out",100,0,plotTMax,100,0,plotTMax);
   h_xIn_deltaX = new TH2F("h_xIn_deltaX",";x_in;x_in-x_out",100,0,plotXMax,100,0,plotXMax);
   h_tIn_deltaT = new TH2F("h_tIn_deltaT",";t_in;t_in-t_out",100,0,plotTMax,100,0,plotTMax);
   h_radMode = new TH1F("h_radMode",";g->gg, g->qq, q->qg",3, -0.5,2.5);
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


 /* if(!CanRadiate(p_in))
    jet.push_back(p_in);
  else{
    Particle p_out1, p_out2;
    Radiate(p_in, p_out1, p_out2);
    DEBUG_MSG("; Particle (type, t, x) = (" 
        << p_in.GetType() << ", " << p_in.GetT() << ", " << p_in.GetX() << ") --> radiated --> ("
        << p_out1.GetType() << ", " << p_out1.GetT() << ", " << p_out1.GetX() << ") and ("
        << p_out2.GetType() << ", " << p_out2.GetT() << ", " << p_out2.GetX() << ")");

    MakeJet(p_out1, jet); // call MakeJet recoursively for daughter particles
    //MakeJet(p_out2, jet);
  }
*/
  Particle p_out1, p_out2;
  if(Radiate(p_in, p_out1, p_out2)){
    DEBUG_MSG("; Particle (type, t, x) = (" 
        << p_in.GetType() << ", " << p_in.GetT() << ", " << p_in.GetX() << ") --> radiated --> ("
        << p_out1.GetType() << ", " << p_out1.GetT() << ", " << p_out1.GetX() << ") and ("
        << p_out2.GetType() << ", " << p_out2.GetT() << ", " << p_out2.GetX() << ")");

    h_xIn_xOut->Fill(p_in.GetX(), p_out1.GetX()); 
    h_tIn_tOut->Fill(p_in.GetT(), p_out1.GetT());
    h_xIn_deltaX->Fill(p_in.GetX(), p_in.GetX()-p_out1.GetX()); 
    h_tIn_deltaT->Fill(p_in.GetT(), p_in.GetT()-p_out1.GetT());

    if(p_in.GetX() == p_out1.GetX()) cout << "x_in == x_out == " << p_in.GetX() << endl;


    MakeJet(p_out1, jet);   // call MakeJet recoursively for daughter particles

    bool followGluons = false;
    if(followGluons)
      MakeJet(p_out2, jet);
    else
      jet.push_back(p_in);
  }
  else
    jet.push_back(p_in);
  return;
}
//-------------------------------------------------
bool FSR::Radiate(Particle p_in, Particle &p_out1, Particle &p_out2)
{
  DEBUG_MSG("FSR::Radiate: " << p_in.GetType());

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
 
    DEBUG_MSG("t_out1_gg " << t_out1_gg << " t_out1_qg " << t_out1_qg );

    if(t_out1_gg < t_out1_qg){ // radiate gluons before quarks
      DEBUG_MSG("\t--> radiating 2 gluons");
      h_radMode->Fill(0);
      t_out1 = t_out1_gg;
      x_out1 = GetXFromP_gg(x_in, r_x);
      type_out1 = gluon;
      type_out2 = gluon;
    }
    else{ // radiate quarks before gluons
      DEBUG_MSG("\t--> radiating 2 quarks");
      h_radMode->Fill(1);
      t_out1 = t_out1_qg;
      x_out1 = GetXFromP_qg(x_in, r_x);
      type_out1 = quark;
      type_out2 = quark;
    }
  }
  else if (p_in.GetType() == quark) { // is a quark
    DEBUG_MSG("quark");
    DEBUG_MSG("\t--> radiating quark and gluon");
    h_radMode->Fill(2);
    
    double r_t_qq = m_rand->Uniform(0,1);
    
    h_rnd->Fill(r_t_qq);
    DEBUG_MSG("r_t_qq " << r_t_qq << ", r_x " << r_x);
    
    t_out1 = GetTFromDelta_qq(t_in, r_t_qq);
    
    x_out1 = GetXFromP_qq(x_in, r_x);
    type_out1 = quark;
    type_out2 = gluon;

    DEBUG_MSG("t_out1 " << t_out1 << " x_out1 " << x_out1 );
  }

  t_out2 = TMath::Sqrt(t_in*t_in - t_out1 *t_out1);
  x_out2 = x_in - x_out1;

  if(t_out1 > m_t0 && t_out2 > m_t0 && x_out1/x_in > m_integrationCutoff && x_out2/x_in > m_integrationCutoff)
    radiated = true;

  DEBUG_MSG("accept radiation");
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
  double C_A = 3.; // FIXME: value
  // numerical cutoff:
  if(z < m_integrationCutoff)
    z = m_integrationCutoff;
  else if(z > (1.-m_integrationCutoff))
    z = 1.-m_integrationCutoff;

  double retVal = C_A * ( z/(1.-z) + (1.-z)/z + z*(1.-z) );
  
  return retVal;
}
//-------------------------------------------------
// Equation (2.114), page 142
double FSR::P_qg(double z)
{
  double T_R = 1./2.; // FIXME: value
  double retVal = T_R * ( z*z + (1.-z)*(1.-z) );

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
  double retVal = C_F * (1.+z*z) / (1.-z);

  return retVal;
}
//-------------------------------------------------
double FSR::Integrate(double (*func)(double), double z0, double z1)
{
  double n = z1/m_precision - z0/ m_precision;
  double integral(0);
  for(int i=0; i<n; i++){
    integral += ( (*func)(z0+m_precision) + (*func)(z0) );
    z0+=m_precision;
  }

  /*if(n<1){
    cout << "z1 " << z1 << " z0 " << z0 << " m_precision " <<  m_precision << " n " << n << endl;
    cout << (0 < n) << endl;
    if(integral == 0) cin.get();
  }*/
  return 0.5 * m_precision * integral;
}

double FSR::IntP_gg(double z0, double z1)
{
//  DEBUG_MSG("FSR::IntP_gg");
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_gg, z0, z1);
//  DEBUG_MSG("\t--> " << retVal);
  return retVal;
}
//-------------------------------------------------
double FSR::IntP_qg(double z0, double z1)
{
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_qg, z0, z1);
  return retVal;
}
//-------------------------------------------------
double FSR::IntP_qq(double z0, double z1)
{
  double retVal = m_alpha_s / TMath::TwoPi() * Integrate(&P_qq, z0, z1);
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
  while(TMath::Exp(-arg) > c && t_out >= m_precision){
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
  while(TMath::Exp(-arg) > c && t_out >= m_precision){
    arg = intP_qg * TMath::Log(t_in/t_out);
    t_out -=  m_precision;
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

  while(TMath::Exp(-arg) > c && t_out >= m_precision){
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
  c = c * IntP_gg(0,1);
  while(integral < c){
    integral += IntP_gg(z0, z0+m_precision);
    z0 += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  double x_out = (z0-1.5*m_precision) * x1;
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
  c = c * IntP_qg(0,1);
  while(integral < c){ 
    integral += IntP_qg(z0, z0+m_precision);
    z0 += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  double x_out = (2*z0-m_precision)/2. * x1;
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
  c = c * IntP_qq(0,1);
    while(integral < c){ 
      integral += IntP_qq(z0, z0 + m_precision);
      z0 += m_precision;
    }
  // upper integration bound is something x0 == x2/x1
  double x_out = (2*z0-m_precision)/2. * x1;
  DEBUG_MSG("\t-->"<<x_out);
  return x_out;
}
//-------------------------------------------------
void FSR::DrawTXPlot(char* pdf){
  TCanvas c("c");
  c.Divide(4,2);
  vector< double > t;
  vector< double > x;
  int n = min((int)m_debugchain.size(),100);
  for(int i=0; i<n; i++){
    t.push_back(m_debugchain[i].GetT());
    x.push_back(m_debugchain[i].GetX());
  }
  TGraph g(n, &t[0], &x[0]);
  if(m_debugchain[0].GetType() == quark)
    g.SetTitle("Initial quark;t;x");
  if(m_debugchain[0].GetType() == gluon)
    g.SetTitle("Initial gluon;t;x");
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

  c.cd(5);
  h_xIn_xOut->Draw("COLZ");
  c.cd(6);
  h_tIn_tOut->Draw("COLZ");
  c.cd(7)->SetLogy();
  h_xIn_deltaX->GetYaxis()->SetRangeUser(1e-9,1);
  h_xIn_deltaX->Draw("COLZ");
  c.cd(8)->SetLogy();
  h_tIn_deltaT->Draw("COLZ");
  c.Print(pdf);


  /*h_t_q->Reset();
  h_x_q->Reset();
  h_t_g->Reset();
  h_x_g->Reset();
  h_rnd->Reset();*/
  h_radMode->Reset();


}


void FSR::DebugPlots(char* pdf, double t_in)
{
  double r_vec[100];
  double t_vec[100];
  double z_vec[100];

  double t_gg_vec[100];
  double t_qg_vec[100];
  double t_qq_vec[100];

  double x_gg_vec[100];
  double x_qg_vec[100];
  double x_qq_vec[100];

  double P_gg_vec[100];
  double P_qg_vec[100];
  double P_qq_vec[100];

  double intP_gg_vec[100];
  double intP_qg_vec[100];
  double intP_qq_vec[100];

  double ratio_intP_gg_vec[100];
  double ratio_intP_qg_vec[100];
  double ratio_intP_qq_vec[100];

  double delta_gg_vec[100];  
  double delta_qg_vec[100];  
  double delta_qq_vec[100];

  double intP_gg = IntP_gg(0, 1);
  double intP_qg = IntP_qg(0, 1);
  double intP_qq = IntP_qq(0, 1);

  double x_in(1);
  for(int i=0; i<100; i++){
    if(i%10 == 0) cout << " loop " << i << endl;
    r_vec[i] = 0.01 + i*0.01;
    t_vec[i] = m_t0 + i*(t_in - m_t0)/100.;
    z_vec[i] = 0.01 + i*0.01;

    t_gg_vec[i] = GetTFromDelta_gg(t_in, r_vec[i]);
    t_qg_vec[i] = GetTFromDelta_qg(t_in, r_vec[i]);
    t_qq_vec[i] = GetTFromDelta_qq(t_in, r_vec[i]);

    x_gg_vec[i] = GetXFromP_gg(x_in, r_vec[i]);
    x_qg_vec[i] = GetXFromP_qg(x_in, r_vec[i]);
    x_qq_vec[i] = GetXFromP_qq(x_in, r_vec[i]);

    P_gg_vec[i] = P_gg(z_vec[i]);
    P_qg_vec[i] = P_qg(z_vec[i]);
    P_qq_vec[i] = P_qq(z_vec[i]);

    intP_gg_vec[i] = IntP_gg(0, z_vec[i]);
    intP_qg_vec[i] = IntP_qg(0, z_vec[i]);
    intP_qq_vec[i] = IntP_qq(0, z_vec[i]);

    ratio_intP_gg_vec[i] = IntP_gg(0, z_vec[i]) / intP_gg;
    ratio_intP_qg_vec[i] = IntP_qg(0, z_vec[i]) / intP_qg;
    ratio_intP_qq_vec[i] = IntP_qq(0, z_vec[i]) / intP_qq;

    delta_gg_vec[i] = Delta_gg(m_t0, t_vec[i]);
    delta_qg_vec[i] = Delta_qg(m_t0, t_vec[i]);
    delta_qq_vec[i] = Delta_qq(m_t0, t_vec[i]);
  }
  TGraph* g_TFromDelta_gg = new TGraph(100,r_vec, t_gg_vec);
  TGraph* g_TFromDelta_qg = new TGraph(100,r_vec, t_qg_vec);
  TGraph* g_TFromDelta_qq = new TGraph(100,r_vec, t_qq_vec);

  TGraph* g_XFromP_gg = new TGraph(100,r_vec, x_gg_vec);
  TGraph* g_XFromP_qg = new TGraph(100,r_vec, x_qg_vec);
  TGraph* g_XFromP_qq = new TGraph(100,r_vec, x_qq_vec);

  TGraph* g_P_gg = new TGraph(100,z_vec, P_gg_vec);
  TGraph* g_P_qg = new TGraph(100,z_vec, P_qg_vec);
  TGraph* g_P_qq = new TGraph(100,z_vec, P_qq_vec);

  TGraph* g_IntP_gg = new TGraph(100,z_vec, intP_gg_vec);
  TGraph* g_IntP_qg = new TGraph(100,z_vec, intP_qg_vec);
  TGraph* g_IntP_qq = new TGraph(100,z_vec, intP_qq_vec);

  TGraph* g_ratio_IntP_gg = new TGraph(100,z_vec, ratio_intP_gg_vec);
  TGraph* g_ratio_IntP_qg = new TGraph(100,z_vec, ratio_intP_qg_vec);
  TGraph* g_ratio_IntP_qq = new TGraph(100,z_vec, ratio_intP_qq_vec);

  TGraph* g_delta_gg = new TGraph(100, t_vec, delta_gg_vec);
  TGraph* g_delta_qg = new TGraph(100, t_vec, delta_qg_vec);
  TGraph* g_delta_qq = new TGraph(100, t_vec, delta_qq_vec);
  
  TCanvas c("c_debug");
  c.Divide(3,2);
  c.cd(1)->SetLogy(1);
  g_TFromDelta_gg->SetLineColor(2);
  g_TFromDelta_gg->Draw("APC");
  g_TFromDelta_gg->GetHistogram()->SetTitle(Form("t_out(r|t_in=%f);r;t_out", t_in));
  g_TFromDelta_gg->GetHistogram()->GetYaxis()->SetRangeUser(1e-1,t_in*1.1);
  g_TFromDelta_qg->SetLineColor(2);
  g_TFromDelta_qg->SetLineStyle(2);
  g_TFromDelta_qg->Draw("PC");
  g_TFromDelta_qq->Draw("PC");
  c.cd(2)->SetLogy(0);
  g_XFromP_gg->SetLineColor(2);
  g_XFromP_gg->Draw("APC");
  g_XFromP_gg->GetHistogram()->SetTitle(Form("x_out(r|x_in=%f);r;x_out", x_in));
  g_XFromP_gg->GetHistogram()->GetYaxis()->SetRangeUser(0,x_in*1.1);
  g_XFromP_qg->SetLineColor(2);
  g_XFromP_qg->SetLineStyle(2);
  g_XFromP_qg->Draw("PC");
  g_XFromP_qq->Draw("PC");
  c.cd(3)->SetLogy(0);
  g_P_qq->Draw("APC");
  g_P_qq->GetHistogram()->SetTitle("P_xx(z);z;P_xx(z)");
  g_P_gg->SetLineColor(2);
  g_P_gg->Draw("PC");
  g_P_qg->SetLineColor(2);
  g_P_qg->SetLineStyle(2);
  g_P_qg->Draw("PC");
  c.cd(4)->SetLogy(1);
  g_IntP_gg->SetLineColor(2);
  g_IntP_gg->Draw("APC");
  g_IntP_gg->GetHistogram()->SetTitle("intP_xx(0,z);z;intP_xx(0,z)");
  g_IntP_gg->GetHistogram()->GetYaxis()->SetRangeUser(1e-4,0.2);
  g_IntP_qg->SetLineColor(2);
  g_IntP_qg->SetLineStyle(2);
  g_IntP_qg->Draw("PC");
  g_IntP_qq->Draw("PC");
  c.cd(5)->SetLogy(0);
  g_ratio_IntP_gg->SetLineColor(2);
  g_ratio_IntP_gg->Draw("APC");
  g_ratio_IntP_gg->GetHistogram()->SetTitle("intP_xx(0,z)/intP_xx(0,1);z;intP_xx(0,z)/intP_xx(0,1)");
  g_ratio_IntP_gg->GetHistogram()->GetYaxis()->SetRangeUser(0,1.1);
  g_ratio_IntP_qg->SetLineColor(2);
  g_ratio_IntP_qg->SetLineStyle(2);
  g_ratio_IntP_qg->Draw("PC");
  g_ratio_IntP_qq->Draw("PC");

  c.cd(5)->SetLogy(0);
  g_delta_gg->SetLineColor(2);
  g_delta_gg->Draw("APC");
  g_delta_gg->GetHistogram()->SetTitle("Delta_xx(m_t0, t);t;Delta_xx(m_t0, t)");
  //  g_delta_gg->GetHistogram()->GetYaxis()->SetRangeUser(0,1.1);
  g_delta_qg->SetLineColor(2);
  g_delta_qq->SetLineStyle(2);
  g_delta_qg->Draw("PC");
  g_delta_qq->Draw("PC");

  c.Print(pdf);
}
