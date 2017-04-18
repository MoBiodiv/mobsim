#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------------------------------------
//simulate spatial ThomaS process
// the function is an efficient re-implementation of the rThomas function from the spatstat package

// [[Rcpp::export]]
DataFrame rThomas_rcpp(int n_points,
                       int n_mother_points,
                       double sigma,
                       double mu,     //mean number of points per cluster,
                       double xmin = 0,
                       double xmax = 1,
                       double ymin = 0,
                       double ymax = 1
)
{
   //simulate mother points
   //double kappa = n_points/mu/((xmax-xmin)*(ymax-ymin)); //density of mother points

   // double expand = 4.0*sigma;
   //
   // double xmin2 = xmin - expand;
   // double xmax2 = xmax + expand;
   //
   // double ymin2 = ymin - expand;
   // double ymax2 = ymax + expand;
   //
   // double lambda_mother = kappa * (xmax2 - xmin2) * (ymax2 - ymin2);

   NumericVector xpoints(n_points);
   NumericVector ypoints(n_points);

   NumericVector xmother;
   NumericVector ymother;

   RNGScope scope;

   //int n_mother_points = as<int>(rpois(1, lambda_mother));
   //int n_mother_points = as<int>(rpois(1, kappa));

   if (n_mother_points > 0){

      // xmother = runif(n_mother_points, xmin2, xmax2);
      // ymother = runif(n_mother_points, ymin2, ymax2);
      xmother = runif(n_mother_points, xmin, xmax);
      ymother = runif(n_mother_points, ymin, ymax);

      double xnew, ynew;
      int imother;

      NumericVector dxy(2);

      for (int ipoint = 0; ipoint < n_points; ++ipoint){

         do {

            imother = as<int>(runif(1, 0, n_mother_points));

            dxy = rnorm(2, 0.0, sigma);
            xnew = xmother[imother] + dxy[0];
            ynew = ymother[imother] + dxy[1];

         } while (!(xnew > xmin && xnew < xmax && ynew > ymin && ynew < ymax));

         xpoints[ipoint] = xnew;
         ypoints[ipoint] = ynew;
      }
   } // end if n_mother_points > 0

   else {
      xpoints = runif(n_points, xmin, xmax);
      ypoints = runif(n_points, ymin, ymax);
   }

   DataFrame xydat = DataFrame::create(_["x"] = xpoints,
                                       _["y"] = ypoints);

   return(xydat);
}
