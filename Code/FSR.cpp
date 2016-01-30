#include "FSR.h"
#include "TMath.h"
#include <iomanip>      // std::setw
#include <algorithm>    // std::min
#include "TGraph.h"
#include "TCanvas.h"

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
  , m_alpha_s(1) //FIXME: Value
  , m_precision(1e-5)
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
  
  if(! CanRadiate(p_in) ) // particle stable, add to final jet
    jet.push_back(p_in);
  else { // particle will continue radiating
    Particle p_out1, p_out2;
    if(Radiate(p_in, p_out1, p_out2)){
      DEBUG_MSG("; Particle (type, t, x) = (" 
          << p_in.GetType() << ", " << p_in.GetT() << ", " << p_in.GetX() << ") --> radiated --> ("
          << p_out1.GetType() << ", " << p_out1.GetT() << ", " << p_out1.GetX() << ") and ("
          << p_out2.GetType() << ", " << p_out2.GetT() << ", " << p_out2.GetX() << ")");

      MakeJet(p_out1, jet); // call MakeJet recoursively for daughter particles
      MakeJet(p_out2, jet);
    }
    else
      cout << "WARNING: Something went wrong p_in didn't radiate while it should have.\n";
  }
  return;
}
//-------------------------------------------------
bool FSR::Radiate(Particle p_in, Particle &p_out1, Particle &p_out2)
{
  DEBUG_MSG("FSR::Radiate");

  double t_out1(0), x_out1(0);
  ParticleType type_out1(undefined), type_out2(undefined);

  double t_in = p_in.GetT();
  double x_in = p_in.GetX();

  double r_x    = m_rand->Uniform(0,1); //TODO: what happens if r_x == 0? What does it mean physically?
  
  DEBUG_MSG("r_x " << r_x);

  bool radiated = false;
  if( p_in.GetType() == gluon ){ // is a gluon
    double r_t_gg = m_rand->Uniform(0,1);
    double r_t_qg = m_rand->Uniform(0,1);
    DEBUG_MSG("r_t_gg " << r_t_gg << " r_t_qg " << r_t_qg);
    if(r_t_gg < r_t_qg && r_t_gg != 0){ // radiate gluons before quarks
      DEBUG_MSG("--> radiating 2 gluons");
      t_out1 = GetTFromDelta_g(m_t0, Delta_g(m_t0, t_in) / r_t_gg);
      x_out1 = GetXFromP_gg(0, x_in, IntP_gg(0,1) * r_x);
      type_out1 = gluon;
      type_out2 = gluon;
      radiated = true;
    }
    else if(r_t_qg != 0){ // radiate quarks before gluons
      DEBUG_MSG("--> radiating 2 quarks");
      t_out1 = GetTFromDelta_g(m_t0, Delta_g(m_t0, t_in) / r_t_qg);
      x_out1 = GetXFromP_qg(0, x_in, IntP_qg(0,1) * r_x);
      type_out1 = quark;
      type_out2 = quark;
      radiated = true;
    }
  }
  else if (p_in.GetType() == quark) { // is a quark
    DEBUG_MSG("--> radiating quark and gluon");
    double r_t_qq = m_rand->Uniform(0,1);
    DEBUG_MSG("r_t_qq " << r_t_qq);
    if(r_t_qq != 0){
      t_out1 = GetTFromDelta_q(m_t0, Delta_q(m_t0, t_in) / r_t_qq);
      x_out1 = GetXFromP_qq(0, x_in, IntP_qq(0,1) * r_x);
      type_out1 = quark;
      type_out2 = gluon;
      radiated = true;
    }
  }

  if(radiated){
    p_out1.SetT(t_out1);
    p_out1.SetX(x_out1);
    p_out1.SetType(type_out1);
    p_out2.SetT(t_in-t_out1); // FIXME: check this
    p_out2.SetX(x_in-x_out1); // FIXME: check this
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
// Equation (2.166), page 157
double FSR::Delta_g(double t0, double t1){
  DEBUG_MSG("FSR::Delta_g");
  double intP_gg = IntP_gg(0,1);
  double intP_qg = IntP_qg(0,1);
  double arg = (intP_gg + intP_qg) * log(t1/t0);
  return exp(-arg);
}
//-------------------------------------------------
// Equation (2.166), page 157
double FSR::Delta_q(double t0, double t1)
{
  DEBUG_MSG("FSR::Delta_q");
  double intP_qq = IntP_qq(0,1);
  return exp( -intP_qq * log(t1/t0) );
}
//-------------------------------------------------
// Equation (2.105), page 140
double FSR::P_gg(double z)
{
  double C_A = 3; // FIXME: value
  // numerical cutoff:
  //if(z<m_precision || z>1-m_precision)
  if(z<0.1 || z>0.9)
    return 0;
  else
    return C_A * ( z/(1-z) + (1-z)/z + z*(1-z) );
}
//-------------------------------------------------
// Equation (2.114), page 142
double FSR::P_qg(double z)
{
  double T_R = 1./2.; // FIXME: value
  return T_R * ( z*z + (1-z)*(1-z) );
}
//-------------------------------------------------
// Equation (2.118), page 143
double FSR::P_qq(double z)
{
  double m_C_F = 4./3.; // FIXME: value
  // numerical cutoff:
  if(z > 1-m_precision)
    return 0;
  else
    return m_C_F * (1+z*z) / (1-z);
}
//-------------------------------------------------
double FSR::IntP_gg(double x0, double x1)
{
  DEBUG_MSG("FSR::IntP_gg");
  // FIXME: integration should definetly done in a smarte way.
  // Ideally also code multiplication should be avoided in IngP_xx
  int n = (x1 - x0) / m_precision;
  double x(0), integral(0);
  for(int i=0; i<n; i++){
    x = x0 + i*m_precision;
    integral += P_gg(x);
  }
  return m_alpha_s / TMath::TwoPi() * m_precision * integral;
}
//-------------------------------------------------
double FSR::IntP_qg(double x0, double x1)
{
  DEBUG_MSG("FSR::IntP_qg");
  // FIXME: integration should definetly done in a smarte way.
  // Ideally also code multiplication should be avoided in IngP_xx
  int n = (x1 - x0) / m_precision;
  double x(0), integral(0);
  for(int i=0; i<n; i++){
    x = x0 + i*m_precision;
    integral += P_qg(x);
  }
  return m_alpha_s / TMath::TwoPi() *  m_precision * integral;
}
//-------------------------------------------------
double FSR::IntP_qq(double x0, double x1)
{
  DEBUG_MSG("FSR::IntP_qq");
  // FIXME: integration should definetly done in a smarte way.
  // Ideally also code multiplication should be avoided in IngP_xx
  int n = (x1 - x0) / m_precision;
  double x(0), integral(0);
  for(int i=0; i<n; i++){
    x = x0 + i*m_precision;
    integral += P_qq(x);
  }
  return m_alpha_s / TMath::TwoPi() * m_precision * integral;
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
double FSR::GetTFromDelta_g(double t_low, double c)
{
  DEBUG_MSG("FSR::GetTFromDelta_g");
  double intP_gg = IntP_gg(0,1);
  double intP_qg = IntP_qg(0,1);

  DEBUG_MSG(intP_gg << " " << intP_qg);

  double arg(0);
  double t0(t_low), t1(t_low);

  // integrate while varying upper boundary until exp(-arg) == c;
  while(exp(-arg) > c){ // FIXME: > or < ?
    arg = (intP_gg + intP_qg) * log(t1/t_low);
    t0 = t1;
    t1 += m_precision;
  }

  return (t0 + t1)/2.;
}
//-------------------------------------------------
// Equation (2.174), page 159
double FSR::GetTFromDelta_q(double t_low, double c){
  double intP_qq = IntP_qq(0,1);

  double arg(0);
  double t0(t_low), t1(t_low);

  // integrate while varying upper boundary until exp(-arg) == c;
  while(exp(-arg) > c){ // FIXME: > or < ?
    arg = intP_qq * log(t1/t_low);
    t0 = t1;
    t1 += m_precision;
  }
  return (t0 + t1)/2.;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_gg(double x0, double x1, double c)
{
  DEBUG_MSG( "FSR::GetXFromR_gg" );
  double integral(0);
  double x(x0);
  c = c/m_precision;
  while(integral < c){
    integral += P_gg(x);
    x0 = x;
    x += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  return x0 * x1;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_qg(double x0, double x1, double c) 
{
  DEBUG_MSG( "FSR::GetXFromR_qg" );
  double integral(0);
  double x(x0);
  c = c/m_precision;
  while(integral < c){ 
    integral += P_qg(x);
    x0 = x;
    x += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  return x0 * x1;
}
//-------------------------------------------------
// Equation (2.175), page 159
double FSR::GetXFromP_qq(double x0, double x1, double c) 
{
  DEBUG_MSG( "FSR::GetXFromR_qg" );
  double integral(0);
  double x(x0);
  c = c/m_precision;
  while(integral < c){
    integral += P_qq(x);
    x0 = x;
    x += m_precision;
  }
  // upper integration bound is something x0 == x2/x1
  return x0 * x1;
}
//-------------------------------------------------
void FSR::DrawTXPlot(char* pdf){
  TCanvas c("c");
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
  g.Draw("A*L");
  m_debugchain.clear();
  c.Print(pdf);
}
