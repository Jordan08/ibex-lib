//============================================================================
//                               I B E X                                   
// File        : ibex_SmearFunction.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : July 19 2012
// Last Update : July 19 2012
//============================================================================

#include "ibex_SmearFunction.h"

using std::pair;

namespace ibex {



  pair<IntervalVector,IntervalVector> SmearFunction::bisect(const IntervalVector& box, int& last_var) {
    IntervalMatrix J(sys.nb_ctr, sys.nb_var);
    int n = box.size();
    sys.f.jacobian(box,J);

    int var = var_to_bisect (J,box);
    if (var == -1)
      return RoundRobin::bisect(box,last_var);
    else
      return box.bisect(var,ratio);
  }

  // computes the variable with the greatest maximal impact
  int SmearMax::var_to_bisect (IntervalMatrix& J, const IntervalVector& box  ) const { 
    double max_magn = NEG_INFINITY;
    int var=-1;
    for (int j=0; j<nbvars; j++) {
      if (box[j].diam()>= w)
	{
	  for (int i=0; i<sys.nb_ctr; i++)
	    
	    if ( J[i][j].mag() * box[j].diam() > max_magn )
	      {
		max_magn = J[i][j].mag()* box[j].diam();
		var = j;
	      }
	}
    }
    return var;
  }


  // computes the variable with the greatest  sum of impacts
  int SmearSum::var_to_bisect(IntervalMatrix& J,const IntervalVector& box  ) const {
    double max_magn = NEG_INFINITY;
    int var = -1;

    for (int j=0; j<nbvars; j++) {
      if (box[j].diam()>= w)
	{
	  double sum_smear=0;
	  for (int i=0; i<sys.nb_ctr; i++) 
	    {
	      sum_smear+= J[i][j].mag() *box[j].diam(); 
	    }
	  if (sum_smear > max_magn)
	    { max_magn = sum_smear;
	      var = j;
	    }
	}
    }
    return var;
  }


  int SmearSumRelative::var_to_bisect(IntervalMatrix & J, const IntervalVector& box ) const {
  double max_magn = NEG_INFINITY;
  int var = -1;
  // the normalizing factor per constraint
  double ctrjsum[sys.nb_ctr];

  for (int i=0; i<sys.nb_ctr; i++) 
    {ctrjsum[i]=0;
      for (int j=0; j<nbvars ; j++)
	{ctrjsum[i]+= J[i][j].mag() * box[j].diam();
	}
    }
      // computes the variable with the maximal sum of normalized impacts
  for (int j=0; j<nbvars; j++) {
    if (box[j].diam()>= w)
      {
	double sum_smear=0;
	for (int i=0; i<sys.nb_ctr; i++) 
	  {if (ctrjsum[i]!=0)
	      sum_smear+= J[i][j].mag() * box[j].diam() / ctrjsum[i]; 
	  }
	if (sum_smear > max_magn)
	  { max_magn = sum_smear;
	    var = j;
	  }
      }
  }
  return var;
}

 



  int SmearMaxRelative::var_to_bisect(IntervalMatrix & J,const IntervalVector& box  ) const {

  double max_magn = NEG_INFINITY;
  int var = -1;

  
  double ctrjsum[sys.nb_ctr]; // the normalizing factor per constraint
  for (int i=0; i<sys.nb_ctr; i++) 
    {ctrjsum[i]=0;
      for (int j=0; j<nbvars ; j++)
	{ctrjsum[i]+= J[i][j].mag() * box[j].diam() ;
	}
    } 

  // computes the variable with the greatest normalized impact
  double maxsmear=0;
  for (int j=0; j<nbvars; j++) {
    if (box[j].diam()>= w)

	for (int i=0; i<sys.nb_ctr; i++) 
	  {
	    if (ctrjsum[i]!=0)
	      maxsmear = J[i][j].mag() * box[j].diam() / ctrjsum[i]; 
	    if (maxsmear > max_magn)
	      { max_magn = maxsmear;
		var = j;
	      }
	  }

  }
  return var;
}



} // end namespace ibex