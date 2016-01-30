/* 
 *  Implementation of QCD Final State Radiation.
 *  This was a rather quick test, but maybe we can follow the general structure
 *
 *  The way I perform integrations is very primitive and doesn't produce good results, in particular for gluon->gluon,gluon
 *  To be honest I was quite surprised I got any results at all
 *
 *  So far I haven't paid any consideration to constants, and variable values and units
 * 
 *  Maybe you can have a look at it (there's a little test in TestFSR.cpp)
 *
 *  Started by Arthur Bolz
 *  Started January 30, 2015
 */
#ifndef __FSR_H__
#define __FSR_H__

#define DEBUG (1)
#if(DEBUG)
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif


#include "iostream"
#include "vector"
#include "TRandom3.h"
#include "TLorentzVector.h"

using namespace std;

enum ParticleType {quark, gluon, undefined};

// I would propose to have our own particle class. Then we don't have to worry about
// what kind of objects other groups work with (TLorentzVectors or whatever) and
// if we need to change something we can easily change it here without having to 
// touch the core FSR code.
class Particle{

  public:

    // simplest constructor
    Particle(ParticleType type=undefined, double t=0, double x=0);
    // if another group works with e.g. TLorentzVector's we can add
    // something like this and then a getter function.
    Particle(ParticleType type, TLorentzVector v);
    ~Particle() {};
    
    void SetType(ParticleType type) { m_type = type; }
    void SetT(double t) { m_t = t;}
    void SetX(double x) { m_x = x;}

    ParticleType GetType(void) { return m_type; }
    double GetT(void) { return m_t; }
    double GetX(void) { return m_x; }


  private:

    ParticleType m_type;
    double m_t;
    double m_x;
};


// final state radiation test class
class FSR{

  public:
    // constructor
    FSR(double t0);
    ~FSR();

    // produce a jet from a particle, call recursively
    void MakeJet(Particle p_in, vector< Particle >& jet);

    void DrawTXPlot(char* pdf);

  private:
    // check whether a particle can still radiate (eg t>t0
    bool CanRadiate(Particle p);
    bool Radiate(Particle p_in, Particle &p_out1, Particle &p_out2);

    double Delta_g(double t0, double t1);
    double Delta_q(double t0, double t1);

    double GetTFromDelta_g(double t_low, double c);
    double GetTFromDelta_q(double t_low, double c);

    double P_gg(double z);
    double P_qg(double z);
    double P_qq(double z);

    double IntP_gg(double x0, double x1);
    double IntP_qg(double x0, double x1);
    double IntP_qq(double x0, double x1);

    double GetXFromP_gg(double x0, double x1, double c);
    double GetXFromP_qg(double x0, double x1, double c);
    double GetXFromP_qq(double x0, double x1, double c);

    TRandom3* m_rand;

    double m_t0;                    // Lower t0 bound, whats a good value, which units do we use (GeV?)
    double m_alpha_s;               // TODO: how to treat alpha_s, at first set constant?, to which value?;
    double m_precision;             // precision of integrations etc. so far integrations are extremely primitive.
    std::vector< Particle > m_debugchain;

};

#endif