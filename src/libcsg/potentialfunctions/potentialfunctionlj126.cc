
#include <votca/csg/potentialfunctions/potentialfunctionlj126.h>

PotentialFunctionLJ126::PotentialFunctionLJ126(const string& name_,const double min_,
                                               const double max_) : PotentialFunction(name_,2,min_,max_){
}

double PotentialFunctionLJ126::CalculateF (const double r) const {

  if ( r >= _min && r <= _cut_off )
    return abs(_lam(0))/pow(r,12) - abs(_lam(1))/pow(r,6) ;
  else
    return 0.0;

}



double PotentialFunctionLJ126::CalculateV (const double r) const {

  if ( r >= _min && r <= _cut_off )
      return 12.0 * abs(_lam(0))/pow(r,12) - 6.0 * abs(_lam(1))/pow(r,6) ;
  else
      return 0.0;
}



// calculate first derivative w.r.t. ith parameter
double PotentialFunctionLJ126::CalculateDF(const int i, const double r) const {

  if ( r >= _min && r <= _cut_off ) {

    switch(i) {
    case 0:
      if ( _lam(0) > 0)
        return 1.0/pow(r,12);
      else
        return -1.0/pow(r,12);
    case 1:
      if ( _lam(1) > 0)
        return -1.0/pow(r,6);
      else 
        return 1.0/pow(r,6);
    }

  }
  return 0.0;
}


// calculate first derivative w.r.t. ith parameter
double PotentialFunctionLJ126::CalculateDV(const int i, const double r) const {

   if ( r >= _min && r <= _cut_off ) {

     switch(i) {
        case 0:
	  if ( _lam(0) > 0)
            return 12.0/pow(r,12);
	  else 
	    return -12.0/pow(r,12);
        case 1:
	  if ( _lam(1) > 0) 
	    return -6.0/pow(r,6);
	  else
	    return 6.0/pow(r,6);
	    
       }
   }
   else
     return 0.0;
}
//

// calculate second derivative w.r.t. ith and jth parameters
double PotentialFunctionLJ126::CalculateD2F(const int i, const int j,
                                            const double r) const {

    return 0.0;

}
