#include <Rcpp.h>
using namespace Rcpp;


double pdfPoisson(unsigned int n, double lambda)
{
  double logSum = std::log(1);
  for(unsigned int i = 2;i <= n; i++)
  {
    logSum += std::log(i);
  }
  double x = -lambda+n*std::log(lambda)-logSum;
  return std::exp(x);
}

double pdfNormal(double x, double mu, double volat)
{
  return std::exp(-(x-mu)*(x-mu)/(2*volat*volat))/(sqrt(2*M_PI)*volat);
}

double pdfMerton(double x, double t, double drift, double volat, double lambda, double a, double b)
{

  double sum = 0;
  for(unsigned int i = 0; i < 200; i++)
  {
    double p = pdfPoisson(i, lambda*t);
    double f = pdfNormal(x, drift*t+i*a, std::sqrt(t*volat*volat+i*b*b));
    sum += p*f;
  }
  return sum;
}


// [[Rcpp::export]]
std::vector<double> pdfMerton(std::vector<double> x, double t, std::vector<double> param)
{
  std::vector<double> sum(x.size());
  double drift = param[0];
  double volat = param[1];
  double lambda = param[2];
  double a = param[3];
  double b = param[4];
  for(unsigned int i = 0; i<x.size(); i++)
  {
    sum[i] = pdfMerton(x[i], t, drift, volat, lambda, a, b);
  }
  return sum;
}
