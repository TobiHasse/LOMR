# include "mex.h"
# include "matrix.h"
# include "math.h"

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void flow_conv_int ( const int N, const double lam_o, const double c[], const double s[], const double ds[], double int_term[] );

/**********************************************************************/

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )

/**********************************************************************/
/*
  Purpose:  

    MEXFUNCTION is a MATLAB/C interface. The original code in this file
    was for computing the Chebyshev U polynomial. 
    The algorithm has been rewritten by Tobias to speed up the flowfield 
    computation in the meander model. Modificatiosn to the file are 
    minimal otherwise and credit for the form of this code goes to the 
    original authors

  Discussion:

    This file should be called "flow_conv_int_trunc.c". (formerly cheby_u.c)
    It should be placed in the MATLAB user's path.  It can either be 
    compiled externally, with a command like

      mex flow_conv_int.c

    creating a compiled MEX file, or, inside of MATLAB, the command

      mex flow_conv_int.c

    accomplishes the same task.  Once the file has been compiled,
    the MATLAB user can invoke the function by typing:

      cx = flow_conv_int ( N, lam_o, C, s, dS )

    The routine computes the flowfield without truncating the weighting 
    function.  The whole upstream curvature series is used.  This can be
    slow computationally for long channels

  Original:
    The original routine calculated the Chebyshev U polynomials and has the
    following citation
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2006
  Author:
    Duane Hanselman and Bruce Littlefield
  Reference:
    Duane Hanselman, Bruce Littlefield,
    Mastering MATLAB 7,
    Pearson Prentice Hall, 2005,
    ISBN: 0-13-143018-1.

    The Mathworks,
    MATLAB External Interfaces,
    The Mathworks, 2000.

  Repurposing
    The format of this code was used by Tobias Ackerman, 2015 as part of 
    dissertation reserach to speed up flowfield computations.

  Parameters:

    Input, int NLHS, the number of output arguments.

    Input, mxArray *PLHS[NLHS], pointers to the output arguments.

    Input, int NRHS, the number of input arguments.

    Input, const mxArray *PRHS[NRHS], pointers to the input arguments.
*/
{
  int    N; 
  double lam_o; 
  double        *c; 
  double        *s; 
  double       *ds; 
  double *int_term;
/*
  GET INPUT:
    From the command line

    int_term = flow_int ( N, lam_o, c[], s[], ds[] )

    Retrieve the first scalar integer input argument N,
    and the second (scalar?) real input argument X.

    cx = cheby_u ( n, x )
*/
  N      = mxGetScalar ( prhs[0] );
  lam_o  = mxGetScalar ( prhs[1] ); 
  c      = mxGetPr     ( prhs[2] );
  s      = mxGetPr     ( prhs[3] );
  ds     = mxGetPr     ( prhs[4] );
/*
  MAKE ROOM FOR OUTPUT:
    Make space for the output argument,
    and copy the pointer to that space.
*/
  plhs[0] = mxCreateDoubleMatrix ( N, 1, mxREAL );

  int_term = mxGetPr ( plhs[0] );
/*
  COMPUTATION:
    Now that we have the interface set up, we can call the C routine.
*/

  flow_conv_int ( N, lam_o, c, s, ds, int_term );

  return;
}
/******************************************************************************/

void flow_conv_int ( const int N, const double lam_o, const double c[], const double s[], const double ds[], double int_term[] )

/******************************************************************************/
/*
  Purpose:  

    Once again, the original purposes of this code have been completely changed
    as I have written an alorithm for the meander model flowfield 
    CHEBY_U evaluates the Chebyshev polynomials of the second kind.

  Original:
    CHEBY_U evaluates the Chebyshev polynomials of the second kind.

  First terms:

    U(0)(X) =   1
    U(1)(X) =   2 X
    U(2)(X) =   4 X^2 -   1
    U(3)(X) =   8 X^3 -   4 X
    U(4)(X) =  16 X^4 -  12 X^2 +  1
    U(5)(X) =  32 X^5 -  32 X^3 +  6 X
    U(6)(X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
    U(7)(X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X

  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2006
  Author:
    John Burkardt
  Parameters:
    Input, int N, the highest polynomial to compute.
    Input, double X, the point at which the polynomials are to be computed.
    Output, double CX[N+1], the values of the N+1 Chebyshev polynomials.

  Flowfield algorithm author:

    Tobias Hasse

  Modified:

    April 17 2015 original source code
    June 3 2021 additional commenting

  MATLAB algorithm, translate into c...
  for ds_idx = 2:N
            up_id = us_idx(ds_idx);
            integrands = C(up_id:ds_idx).*exp(lam_o*(s(ds_idx)-s(up_id:ds_idx)));
            int_term(ds_idx) = dS(up_id:ds_idx-1)' * (integrands(1:end-1) + integrands(2:end))/2;
  end
 
 */
{
  
  int i;
  int j;
  double ig;
  double ig_old;
/* loop for testing inputs */
//   for ( i = 0; i <N; i++ )
//   {
//       mexPrintf("Hello World! %u \n", us_idx[i]);
//       int_term[i] = i;
//   }

// mexPrintf("hello World\n");
/* not truncated version */
for ( i = 1; i <N; i++ )
  {
      ig = c[0] * exp( lam_o * ( s[i]-s[0] )); 
      for ( j = 1; j <= i; j++ )
      {
          ig_old = ig;
          ig = c[j] * exp( lam_o * ( s[i]-s[j] ) );
          int_term[i] = int_term[i] + ds[j-1] * ( ig_old + ig ) / 2;
      }
  }
  return;
}
 