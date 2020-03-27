#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------------------------------------
   //simulate spatial ThomaS process
// the function is an efficient re-implementation of the rThomas function from the spatstat package

// [[Rcpp::export]]
DataFrame rThomas_rcpp(int n_points,
                       int n_mother_points,
                       NumericVector xmother,
                       NumericVector ymother,
                       double sigma,
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
   
   
   RNGScope scope;
   
   //int n_mother_points = as<int>(rpois(1, lambda_mother));
   //int n_mother_points = as<int>(rpois(1, kappa));
   
	bool mother_points_specified = Rcpp::na_omit(xmother).size() > 0;
	
   if (n_mother_points > 0 || mother_points_specified) {
      if(!LogicalVector::is_na(n_mother_points) & all(is_na(xmother)) & all(is_na(ymother))) {
			// xmother = runif(n_mother_points, xmin2, xmax2);
			// ymother = runif(n_mother_points, ymin2, ymax2);
			xmother = runif(n_mother_points, xmin, xmax);
			ymother = runif(n_mother_points, ymin, ymax);
      } else {
			n_mother_points = xmother.size();
		}
		
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
   } // end if n_mother_points > 0 || xmother.size() > 0 : if clustering
   
   else {
      xpoints = runif(n_points, xmin, xmax);
      ypoints = runif(n_points, ymin, ymax);
   }
   
   DataFrame xydat = DataFrame::create(_["x"] = xpoints,
                                        _["y"] = ypoints);
   
   return(xydat);
}
