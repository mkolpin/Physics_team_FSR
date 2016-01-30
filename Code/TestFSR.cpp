#include "FSR.h"

#include "iostream"
#include "math.h"
#include <iomanip>      // std::setw
#include <TH1.h>

using namespace std;

double Integral( double(*f)(double), double x_low, double x_high, double delta_x=1e-6){
  int n = (x_high - x_low) / delta_x;
  double x(0), integral(0);
  for(int i=0; i<n; i++){
    x = x_low + i*delta_x;
    integral += (*f)(x);
  }
  return delta_x * integral;
}

double GetXHigh( double(*f)(double), double x_low, double C, double delta_x=1e-6){
  double x_0(x_low), x_1(x_low);
  double integral=0;
  C = C/delta_x;
  while(integral < C){
    integral += (*f)(x_1);
    x_0 = x_1;
    x_1 += delta_x;
  }

  return (x_0 + x_1)/2.;  
}

double testFCN(double x)
{
  return x;
}
double integralTestFCN(double x_low, double x_high)
{
  return 0.5*(x_high*x_high - x_low*x_low); 
}

double xHighTestFCN(double x_low, double C)
{
  return sqrt(2*C + x_low*x_low);
}

int main()
{
  double x_low(0);
  double x_high(0);
  double integral_true(0), integral_est(0);
  cout << "x_high          integral_true   integral_est    x_high_est     \n";
  for(int i=1; i<=10; i++){
    x_high = i*0.5;
    integral_est = Integral(testFCN, x_low, x_high);
    integral_true = integralTestFCN(x_low, x_high);
    cout << setw(10) << x_high << " " 
      << setw(10) << integral_true << " "
      << setw(10) << integral_est << " "
      << setw(10) << GetXHigh(testFCN, x_low, integral_true) << endl;
  }


  // FSR Test
  FSR test(1e-2);
  Particle p_quark(quark, .9, .9);
  Particle p_gluon(gluon, .9, .9);
  vector< Particle > jet_out;
  test.MakeJet(p_quark, jet_out);
  test.DrawTXPlot("TestFSR.pdf(");
  jet_out.clear();
  test.MakeJet(p_gluon, jet_out);
  test.DrawTXPlot("TestFSR.pdf)");


  return 1;
}
