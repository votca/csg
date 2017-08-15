
#include <votca/csg/potentialfunctions/potentialfunctioncbspl.h>

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const string& name_,const int nlam_,
                                               const double min_, const double max_) :
  PotentialFunction(name_,nlam_,min_,max_) {

  /* Here nlam_ is the total number of coeff values that are to be optimized
   * To ensure that potential and force go to zero smoothly near cut-off,
   * as suggested in Ref. PCCP, 11, 1901, 2009, coeff values leading up to
   * cut-off and beyond take a value of zero.
   *
   * Since region less than rmin is not sampled sufficiently for stability
   * first _nexcl coefficients are not optimized instead their values are
   * extrapolated from first statistically significant knot values near rmin
   */

  int nknots;

  nknots = _lam.size();

  _nbreak = nknots - 2;

  _dr = ( _cut_off )/( double (_nbreak - 1) );

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak.resize(nknots,false);
  _rbreak.clear();

  for( int i = 0; i < nknots; i++)
    _rbreak(i) =  i*_dr;

  // exclude knots corresponding to r <= _min
  _nexcl = min( int( ( _min )/_dr ), _nbreak - 2 ) + 1;

  // account for finite numerical division of _min/_dr
  // e.g. 0.24/0.02 may result in 11.99999999999999
  if( _rbreak(_nexcl) == _min ) _nexcl++;

  // fixing last 4 knots to zeros is reasonable
  _ncutcoeff = 4;

  // check if we have enough parameters to optimize
  if((int(_lam.size()) - _nexcl - _ncutcoeff) < 1)
    {
      throw std::runtime_error("In potential "+_name+": no parameters to optimize!\n"
                               "All the knot values fall in the range of either excluded (due to high repulsive region) or cut-off region.\n"
                               "This issue can be resolved by one or combination of following steps:\n"
                               "1. Make sure you are using large-enough cut-off for this CG potential.\n"
                               "2. Make sure the CG-MD runs are sufficiently long and CG-MD RDF are statistically reliable.\n"
                               "3. Use more knot values.\n");
    }

  _M.resize(4,4,false);
  _M.clear();
  _M(0,0) =  1.0; _M(0,1) =  4.0; _M(0,2) =  1.0; _M(0,3) = 0.0;
  _M(1,0) = -3.0; _M(1,1) =  0.0; _M(1,2) =  3.0; _M(1,3) = 0.0;
  _M(2,0) =  3.0; _M(2,1) = -6.0; _M(2,2) =  3.0; _M(2,3) = 0.0;
  _M(3,0) = -1.0; _M(3,1) =  3.0; _M(3,2) = -3.0; _M(3,3) = 1.0;
  _M /= 6.0;

}

int PotentialFunctionCBSPL::getOptParamSize() const {

  return int(_lam.size()) - _nexcl - _ncutcoeff;

}

void PotentialFunctionCBSPL::setParam(string filename) {

  Table param;
  param.Load(filename);
  _lam.clear();

  if( param.size() != _lam.size()) {

    throw std::runtime_error("In potential "+_name+": parameters size mismatch!\n"
                             "Check input parameter file \""
                             + filename + "\" \nThere should be "
                             + boost::lexical_cast<string>( _lam.size() ) + " parameters");
  } else {

    // force last _ncutcoeff to zero
    for( unsigned int i = 0; i < _lam.size() - _ncutcoeff; i++)
      _lam(i) = param.y(i);

  }

}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  extrapolExclParam();

  Table param;
  param.SetHasYErr(false);
  param.resize(_lam.size(), false);

  // write extrapolated knots with flag 'o'
  // points close to rmin can also be stastically not reliable
  // so flag 3 more points next to rmin as 'o'
  for (int i = 0; i < _nexcl+3; i++)
    param.set(i, _rbreak(i), _lam(i), 'o');

  for (unsigned int i = _nexcl+3; i < _lam.size() ; i++)
    param.set(i, _rbreak(i), _lam(i), 'i');

  param.Save(filename);

}

void PotentialFunctionCBSPL::SavePotTab(const string& filename,
                                        const double step,
                                        const double rmin,
                                        const double rcut)
{
  extrapolExclParam();
  PotentialFunction::SavePotTab(filename,step,rmin,rcut);
}

void PotentialFunctionCBSPL::SavePotTab(const string& filename,
                                        const double step)
{
  extrapolExclParam();
  PotentialFunction::SavePotTab(filename,step);
}

void PotentialFunctionCBSPL::extrapolExclParam(){

  double u0 = _lam(_nexcl);
  double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
    (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
  double r0 = _rbreak(_nexcl);

  /* If the slope m is positive then the potential core
   * will be attractive. So, artificially forcing core to be
   * repulsive by setting m = -m
   */
  if( m > 0)
    {
      cout << _name << " potential's extrapolated core is attractive!" << endl;
      cout << "Artifically enforcing repulsive core.\n" << endl;
      m *= -1.0;
    }
  // using linear extrapolation
  // u(r) = ar + b
  // a = m
  // b = - m*r0 + u0
  // m = (u1-u0)/(r1-r0)

  double a = m;
  double b = -1.0*m*r0 + u0;
  for (int i = 0; i < _nexcl; i++)
    _lam(i) = a*_rbreak(i) + b;

}

void PotentialFunctionCBSPL::setOptParam(const int i, const double val){

  _lam( i + _nexcl ) = val;

}

double PotentialFunctionCBSPL::getOptParam(const int i) const{

  return _lam( i + _nexcl );

}

double PotentialFunctionCBSPL::CalculateF (const double r) const {

  if( r <= _cut_off){

    double u = 0.0;

    ub::vector<double> R;
    ub::vector<double> B;
    R.resize(4,false); R.clear();
    B.resize(4,false); B.clear();

    int indx = min( int( r /_dr ) , _nbreak-2 );;
    double rk = indx*_dr;;
    double t = (r - rk)/_dr;

    R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;
    ub::vector<double> RM = ub::prod(R,_M);
    B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
    B(3) = _lam(indx+3);

    u += ub::inner_prod(B,RM);

    return u;

  } else
    return 0.0;

}
double PotentialFunctionCBSPL::CalculateV (const double r) const {

   return 0.0;
}
// calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const int i, const double r) const{

  // since first _nexcl parameters are not optimized for stability reasons


  if ( r <= _cut_off ) {

    unsigned int i_opt = i + _nexcl;
    unsigned int indx;
    double rk;

    indx = min( int( ( r )/_dr ), _nbreak-2 );
    rk = indx*_dr;

    if ( i_opt >= indx && i_opt <= indx+3 ){

      ub::vector<double> R;
      R.resize(4,false); R.clear();

      double t = ( r - rk)/_dr;

      R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

      ub::vector<double> RM = ub::prod(R,_M);

      return RM(i_opt-indx);

    }else
      return 0.0;

  } else
    return 0.0;

}
double PotentialFunctionCBSPL::CalculateDV(const int i, const double r) const{
 if ( r <= _cut_off ) {

     unsigned int i_opt = i + _nexcl;
     unsigned int indx;
     double rk;
     indx = min( int( ( r )/_dr ), _nbreak-2 );
     rk = indx*_dr;
     
     if ( i_opt >= indx && i_opt <= indx+3 ){
     
     ub::vector<double> R;
     R.resize(4,false); R.clear();

     double t = ( r - rk)/_dr;

     R(0) = 0.0; R(1) = - 1.0*r/_dr; R(2) = - 2.0*(r - rk)*r/(_dr*_dr); R(3) = -3.0*(r-rk)*(r-rk)*r/(_dr*_dr*_dr);

     ub::vector<double> RM = ub::prod(R,_M);
     
     return RM(i_opt-indx);
     
     }else
     
     return 0.0;
     
     } else 
 return 0.0;
}
// calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const int i, const int j,
                                            const double r) const {

  return 0.0;

}
