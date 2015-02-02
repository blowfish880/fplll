#include <fplll.h>
#include <chrono>
#include <math.h> 

using namespace std;
using namespace fplll;

double timedReduce(IntMatrix& M, BKZParam& param) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  IntMatrix U;
  int re = bkzReduction(&M, &U, param);
  
  end = std::chrono::system_clock::now();
  
  if ( re != RED_SUCCESS ) {
    cerr << getRedStatusStr(re) << endl;
    abort();
  }
  return std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;
}

int main(int argc, char *argv[]) {
  // first argument is the seed for the random lattice
  RandGen::initWithSeed(atoi(argv[1]));
  // second argument is the lattice dimension
  int m = atoi(argv[2]);
  
  // generate random hard lattice
  ZZ_mat<mpz_t> M(m, m+1);
  M.gen_intrel(m*100);
  
  // and LLL reduce it
  lllReduction(M);
  
  // parameters for BKZ - beta=m for HKZ
  int beta = m;
  BKZParam param;
  param.blockSize = beta;
  param.flags = BKZ_DEFAULT;
  
  // uncomment for linear pruning
  //~ param.enableLinearPruning(beta);
  
  // recursive kannan style preprocessing
  
  BKZParam* tmp = &param;
  beta = m;
  while (beta > 30) {
      tmp->kan = new BKZParam;
      tmp = tmp->kan;
      // block size for recursive call
      beta -= ilogb(beta);
      tmp->blockSize = beta;
      // uncomment for linear pruning in recursive call
      //~ tmp->enableLinearPruning(beta);
  }
  tmp->kan = NULL;
  
  
  cout << timedReduce(M, param) << endl;
  
  return 0;
}
