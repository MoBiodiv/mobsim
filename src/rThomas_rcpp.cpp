#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------------------------------------
//simulate spatial ThomaS process
// the function is an efficient re-implementation of the rThomas function from the spatstat package

// [[Rcpp::export]]
DataFrame rThomas_rcpp(int nPoints,
                       double sigma,
                       double mu,     //mean number of points per cluster,
                       double xmin = 0,
                       double xmax = 1,
                       double ymin = 0,
                       double ymax = 1
)
{
   //simulate mother points
   double kappa = nPoints/mu/((xmax-xmin)*(ymax-ymin)); //density of mother points

   // double expand = 4.0*sigma;
   //
   // double xmin2 = xmin - expand;
   // double xmax2 = xmax + expand;
   //
   // double ymin2 = ymin - expand;
   // double ymax2 = ymax + expand;
   //
   // double lambda_mother = kappa * (xmax2 - xmin2) * (ymax2 - ymin2);

   NumericVector xpoints(nPoints);
   NumericVector ypoints(nPoints);

   NumericVector xmother;
   NumericVector ymother;

   RNGScope scope;

   //int nMotherPoints = as<int>(rpois(1, lambda_mother));
   int nMotherPoints = as<int>(rpois(1, kappa));

   if (nMotherPoints > 0){

      // xmother = runif(nMotherPoints, xmin2, xmax2);
      // ymother = runif(nMotherPoints, ymin2, ymax2);
      xmother = runif(nMotherPoints, xmin, xmax);
      ymother = runif(nMotherPoints, ymin, ymax);

      double xnew, ynew;
      int imother;

      NumericVector dxy(2);

      for (int ipoint = 0; ipoint < nPoints; ++ipoint){

         do {

            imother = as<int>(runif(1, 0, nMotherPoints));

            dxy = rnorm(2, 0.0, sigma);
            xnew = xmother[imother] + dxy[0];
            ynew = ymother[imother] + dxy[1];

         } while (!(xnew > xmin && xnew < xmax && ynew > ymin && ynew < ymax));

         xpoints[ipoint] = xnew;
         ypoints[ipoint] = ynew;
      }
   } // end if nMotherPoints > 0

   else {
      xpoints = runif(nPoints, xmin, xmax);
      ypoints = runif(nPoints, ymin, ymax);
   }

   DataFrame xydat = DataFrame::create(_["x"] = xpoints,
                                       _["y"] = ypoints);

   return(xydat);
}



