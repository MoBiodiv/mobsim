#include <Rcpp.h>
using namespace Rcpp;

//' Thomas process individual distribution simulation for one species
//'
//' Usually used internally inside \code{\link{sim_thomas_coords}}
//' This function randomly draws points (individuals) around one or several mother points using Rcpp.
//' The function is an efficient re-implementation of the rThomas function from the spatstat package.
//'
//' @name rThomas_rcpp
//'
//' @param n_points The total number of points (individuals).
//' @param n_mother_points Number of mother points (= cluster centres).
//' @param xmother Vector of \code{n_mother_points} x coordinates for the mother points.
//' @param ymother Vector of \code{n_mother_points} y coordinates for the mother points.
//' @param sigma Mean displacement (along each coordinate axes) of a point from
//' its mother point (= cluster centre).
//' @param xmin Left limit, \code{default}=0.
//' @param xmax Right limit, \code{default}=1.
//' @param ymin Bottom limit, \code{default}=0.
//' @param ymax Top limit, \code{default}=1.
//'
//' @return A dataframe with x and y coordinates.
//'
//' @author Felix May, Alban Sagouis
//' @export


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

   NumericVector xpoints(n_points);
   NumericVector ypoints(n_points);

   RNGScope scope;

	bool mother_points_specified = Rcpp::na_omit(xmother).size() > 0;

   // Original version
   // if (n_mother_points > 0 || mother_points_specified) {	// if clustering
   //   if (!LogicalVector::is_na(n_mother_points) &&
   //       all(is_na(xmother)) & all(is_na(ymother))) {

   bool is_na_x_mother = all(is_na(xmother)).is_true();
   bool is_na_y_mother = all(is_na(ymother)).is_true();

   if (n_mother_points > 0 || mother_points_specified) {	// if clustering
      // n_mother_points is an integer not a logical.
      // When sim_thomas_coords is coded correctly, it should never be NA I guess.
      // So I remove it here
      if (is_na_x_mother && is_na_y_mother) {
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
   } // end if n_mother_points > 0 || xmother.size() > 0

   else {
      xpoints = runif(n_points, xmin, xmax);
      ypoints = runif(n_points, ymin, ymax);
   }

   DataFrame xydat = DataFrame::create(_["x"] = xpoints,
                                       _["y"] = ypoints);

   return(xydat);
}
