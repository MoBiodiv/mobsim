#include <algorithm>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;

class CNeighbourDist
{
public:
   double distance;
   //String  ID_spec;
   int ID_spec;

   CNeighbourDist() : distance(0), ID_spec(0) {}

   bool operator< (const CNeighbourDist& neighb1) const
   {
      return (distance < neighb1.distance);
   }
};


//distance function
double DistXY(double x1, double y1, double x2, double y2)
{
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

// [[Rcpp::export]]
NumericVector sSAC1_C(NumericVector x,
                      NumericVector y,
                      IntegerVector id_spec
                      )
{
   int N = x.size();
   NumericVector sSAC(N);

   CNeighbourDist* NeighDistVec = new CNeighbourDist[N];

   std::map<int,int> SpecAbund;

   for (int i=0; i < N; i++){
      // get distance for each neighbour
      for (int j=0; j < N; j++){
         NeighDistVec[j].distance = DistXY(x[i],y[i],x[j],y[j]);
         NeighDistVec[j].ID_spec  = id_spec[j];
      }
      std::sort(NeighDistVec,NeighDistVec+N);

      SpecAbund.clear();
      for (int k=0; k < N; k++){
         ++SpecAbund[NeighDistVec[k].ID_spec];
         sSAC[k] = sSAC[k] + SpecAbund.size();
      }

      Rcpp::checkUserInterrupt();
   }

   for (int i=0; i<N; i++)
      sSAC[i] = (double) sSAC[i]/N;

   return(sSAC);
}
