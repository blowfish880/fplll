/* Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2011 Xavier Pujol.

   This file is part of fplll. fplll is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   fplll is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with fplll. If not, see <http://www.gnu.org/licenses/>. */

#include "wrapper.h"
#include "fplll.h"

#include "gso.cpp"
#include "lll.cpp"
#include "wrapper.cpp"
#include "enumerate.cpp"
#include "svpcvp.cpp"
#include "bkz.cpp"

FPLLL_BEGIN_NAMESPACE

template<class ZT>
void zerosFirst(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInvT) {
  int i, d = b.getRows();
  for (i = d; i > 0 && b[i - 1].is_zero(); i--) {}
  if (i > 0 && i < d) {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!uInvT.empty())
      uInvT.rotate(0, i, d - 1);
  }
}

template<class ZT>
void zerosLast(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInvT) {
  int i, d = b.getRows();
  for (i = 0; i < d && b[i].is_zero(); i++) {}
  if (i > 0 && i < d) {
    b.rotate(0, i, d - 1);
    if (!u.empty())
      u.rotate(0, i, d - 1);
    if (!uInvT.empty())
      uInvT.rotate(0, i, d - 1);
  }
}

template<class ZT, class FT>
int lllReductionZF(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
                   double delta, double eta, LLLMethod method, int flags) {
  int gsoFlags = 0;
  if (b.getRows() == 0 || b.getCols() == 0) return RED_SUCCESS;
  if (method == LM_PROVED) gsoFlags |= GSO_INT_GRAM;
  if (method == LM_FAST)   gsoFlags |= GSO_ROW_EXPO | GSO_OP_FORCE_LONG;
  MatGSO<Z_NR<ZT>, FP_NR<FT> > mGSO(b, u, uInv, gsoFlags);
  LLLReduction<Z_NR<ZT>, FP_NR<FT> > lllObj(mGSO, delta, eta, flags);
  lllObj.lll();
  return lllObj.status;
}

template<class ZT>
int lllReductionWrapper(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
        double delta, double eta, FloatType floatType,
        int precision, int flags) {
  FPLLL_ABORT("The wrapper method works only with integer type mpz");
  return RED_LLL_FAILURE;
}

template<>
int lllReductionWrapper(IntMatrix& b, IntMatrix& u, IntMatrix& uInv,
        double delta, double eta, FloatType floatType,
        int precision, int flags) {
  FPLLL_CHECK(floatType == FT_DEFAULT,
      "The floating point type cannot be specified with the wrapper method");
  FPLLL_CHECK(precision == 0,
      "The precision cannot be specified with the wrapper method");
  Wrapper wrapper(b, u, uInv, delta, eta, flags);
  wrapper.lll();
  zerosFirst(b, u, uInv);
  return wrapper.status;
}

template<class ZT>
int lllReductionZ(ZZ_mat<ZT>& b, ZZ_mat<ZT>& u, ZZ_mat<ZT>& uInv,
                  double delta, double eta,
                  LLLMethod method, IntType intType, FloatType floatType,
                  int precision, int flags) {

  if (method == LM_WRAPPER)
    return lllReductionWrapper(b, u, uInv, delta, eta,
            floatType, precision, flags);

  FPLLL_CHECK(!(method == LM_PROVED && (flags & LLL_EARLY_RED)),
    "LLL method 'proved' with early reduction is not implemented");

  // Computes the parameters required for the proved version
  int goodPrec = l2MinPrec(b.getRows(), delta, eta, LLL_DEF_EPSILON);

  // Sets the parameters and checks the consistency
  int selPrec = 0;
  if (method == LM_PROVED) {
    selPrec = (precision != 0) ? precision : goodPrec;
  }
  else {
    selPrec = (precision != 0) ? precision : PREC_DOUBLE;
  }

  FloatType selFT = floatType;

  if (precision != 0) {
    if (selFT == FT_DEFAULT) {
      selFT = FT_MPFR;
    }
    FPLLL_CHECK(selFT == FT_MPFR,
            "The floating type must be mpfr when the precision is specified");
  }

  if (selFT == FT_DEFAULT) {
    if (method == LM_FAST)
      selFT = FT_DOUBLE;
#ifdef FPLLL_WITH_DPE
    else if (selPrec <= static_cast<int>(FP_NR<dpe_t>::getprec()))
      selFT = FT_DPE;
#endif
    else
      selFT = FT_MPFR;
  }
  else if (method == LM_FAST
           && (selFT != FT_DOUBLE && selFT != FT_LONG_DOUBLE)) {
    FPLLL_ABORT("'double' or 'long double' required for "
                  << LLL_METHOD_STR[method]);
  }

  if (selFT == FT_DOUBLE)
    selPrec = FP_NR<double>::getprec();
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE)
    selPrec = FP_NR<long double>::getprec();
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE)
    selPrec = FP_NR<dpe_t>::getprec();
#endif

  if (flags & LLL_VERBOSE) {
    cerr << "Starting LLL method '" << LLL_METHOD_STR[method] << "'" << endl
         << "  integer type '" << INT_TYPE_STR[intType] << "'" << endl
         << "  floating point type '" << FLOAT_TYPE_STR[selFT] << "'" << endl;
    if (method != LM_PROVED || intType != ZT_MPZ || selFT == FT_DOUBLE)
      cerr << "  The reduction is not guaranteed";
    else if (selPrec < goodPrec)
      cerr << "  prec < " << goodPrec << ", the reduction is not guaranteed";
    else {
      cerr << "  prec >= " << goodPrec << ", the reduction is guaranteed";
    }
    cerr << endl;
  }

  // Applies the selected method
  int status;
  if (selFT == FT_DOUBLE) {
    status = lllReductionZF<ZT, double>(b, u, uInv,
            delta, eta, method, flags);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE) {
    status = lllReductionZF<ZT, long double>(b, u, uInv,
            delta, eta, method, flags);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE) {
    status = lllReductionZF<ZT, dpe_t>(b, u, uInv,
            delta, eta, method, flags);
  }
#endif
  else if (selFT == FT_MPFR) {
    int oldPrec = FP_NR<mpfr_t>::setprec(selPrec);
    status = lllReductionZF<ZT, mpfr_t>(b, u, uInv,
            delta, eta, method, flags);
    FP_NR<mpfr_t>::setprec(oldPrec);
  }
  else {
    FPLLL_ABORT("Compiled without support for LLL reduction with " 
                 << FLOAT_TYPE_STR[selFT]);
  }
  zerosFirst(b, u, uInv);
  return status;
}

/* We define LLL for each input type instead of using a template,
   in order to force the compiler to instantiate the functions. */

#define FPLLL_DEFINE_LLL(T, idT)                                        \
int lllReduction(ZZ_mat<T>& b, double delta, double eta,                \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  ZZ_mat<T> emptyMat; /* Empty u -> transform disabled */               \
  return lllReductionZ<T>(b, emptyMat, emptyMat, delta, eta, method,    \
            idT, floatType, precision, flags);                          \
}                                                                       \
                                                                        \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u, double delta, double eta,  \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  ZZ_mat<T> emptyMat;                                                   \
  if (u.empty()) u.gen_identity(b.getRows());                           \
  return lllReductionZ<T>(b, u, emptyMat, delta, eta, method,           \
            idT, floatType, precision, flags);                          \
}                                                                       \
                                                                        \
int lllReduction(ZZ_mat<T>& b, ZZ_mat<T>& u, ZZ_mat<T>& uInv,           \
                 double delta, double eta,                              \
                 LLLMethod method, FloatType floatType,                 \
                 int precision, int flags) {                            \
  if (u.empty()) u.gen_identity(b.getRows());                           \
  if (uInv.empty()) uInv.gen_identity(b.getRows());                     \
  uInv.transpose();                                                     \
  int status = lllReductionZ<T>(b, u, uInv, delta, eta, method,         \
            idT, floatType, precision, flags);                          \
  uInv.transpose();                                                     \
  return status;                                                        \
}

FPLLL_DEFINE_LLL(mpz_t, ZT_MPZ)

#ifdef FPLLL_WITH_ZLONG
FPLLL_DEFINE_LLL(long, ZT_LONG)
#endif

#ifdef FPLLL_WITH_ZDOUBLE
FPLLL_DEFINE_LLL(double, ZT_DOUBLE)
#endif

template<class FT>
int bkzReductionF(IntMatrix &b, const BKZParam& param, int selFT, double lllDelta, IntMatrix& u, IntMatrix& uInv) {
  int gsoFlags = 0;
  if (b.getRows() == 0 || b.getCols() == 0)
    return RED_SUCCESS;
  if (selFT == FT_DOUBLE || selFT == FT_LONG_DOUBLE)
    gsoFlags |= GSO_ROW_EXPO;
  MatGSO<Integer, FT> mGSO(b, u, uInv, gsoFlags);
  LLLReduction<Integer, FT> lllObj(mGSO, lllDelta, LLL_DEF_ETA, LLL_DEFAULT);
  BKZReduction<FT> bkzObj(mGSO, lllObj, param);
  bkzObj.bkz();
  return bkzObj.status;
}

int bkzReduction(IntMatrix* B, IntMatrix *U, const BKZParam& param, FloatType floatType, int precision) {
  IntMatrix emptyMat;
  IntMatrix& u = U ? *U : emptyMat;
  IntMatrix& uInv = emptyMat;
  FPLLL_CHECK(B, "B == NULL in bkzReduction");

  if (U && u.empty()) {
    u.gen_identity(B->getRows());
  }

  double lllDelta = param.delta < 1 ? param.delta : LLL_DEF_DELTA;

  FloatType selFT = (floatType != FT_DEFAULT) ? floatType : FT_DOUBLE;
  FPLLL_CHECK(!(selFT == FT_MPFR && precision == 0),
               "Missing precision for BKZ with floating point type mpfr");

  if (param.flags & BKZ_NO_LLL)
    zerosLast(*B, u, uInv);
  else {
    Wrapper wrapper(*B, u, uInv, lllDelta, LLL_DEF_ETA, LLL_DEFAULT);
    if (!wrapper.lll()) return wrapper.status;
  }

  int status;
  if (selFT == FT_DOUBLE) {
    status = bkzReductionF< FP_NR<double> >(*B, param, selFT, lllDelta, u, uInv);
  }
#ifdef FPLLL_WITH_LONG_DOUBLE
  else if (selFT == FT_LONG_DOUBLE) {
    status = bkzReductionF< FP_NR<long double> >(*B, param, selFT, lllDelta, u, uInv);
  }
#endif
#ifdef FPLLL_WITH_DPE
  else if (selFT == FT_DPE) {
    status = bkzReductionF< FP_NR<dpe_t> >(*B, param, selFT, lllDelta, u, uInv);
  }
#endif
  else if (selFT == FT_MPFR) {
    int oldPrec = FP_NR<mpfr_t>::setprec(precision);
    status = bkzReductionF< FP_NR<mpfr_t> >(*B, param, selFT, lllDelta, u, uInv);
    FP_NR<mpfr_t>::setprec(oldPrec);
  }
  else {
    FPLLL_ABORT("Compiled without support for BKZ reduction with " 
                 << FLOAT_TYPE_STR[selFT]);
  }
  zerosFirst(*B, u, uInv);
  return status;
}

int bkzReduction(IntMatrix& b, int blockSize, int flags, FloatType floatType, int precision) {
  BKZParam param;
  param.blockSize = blockSize;
  param.flags = flags;
  return bkzReduction(&b, NULL, param, floatType, precision);
}

int bkzReduction(IntMatrix& b, IntMatrix& u, int blockSize, int flags, FloatType floatType, int precision) {
  BKZParam param;
  param.blockSize = blockSize;
  param.flags = flags;
  return bkzReduction(&b, &u, param, floatType, precision);
}

int hkzReduction(IntMatrix& b, int flags) {
  BKZParam param;
  param.blockSize = b.getRows();
  param.delta = 1;
  if (flags & HKZ_VERBOSE) param.flags |= BKZ_VERBOSE;
  return bkzReduction(&b, NULL, param);
}

const char* getRedStatusStr(int status) {
  if (status >= 0 && status < RED_STATUS_MAX)
    return RED_STATUS_STR[status];
  else
    return "unknown error";
}


int dSvpReduce(IntMatrix& b, int start, int end) {
  int gsoFlags = 0;
  if (b.getRows() == 0 || b.getCols() == 0)
    return RED_SUCCESS;
  IntMatrix emptyMat;
  //~ cout << "GSO object" << endl;
  MatGSO<Integer, Float> mGSO(b, emptyMat, emptyMat, gsoFlags);
  //~ cout << "LLL" << endl;
  LLLReduction<Integer, Float> lllObj(mGSO, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
  BKZParam param;
  param.blockSize = end - start;
  param.delta = LLL_DEF_DELTA;
  
  //~ cout << "BKZ object" << endl;
  BKZReduction<Float> bkzObj(mGSO, lllObj, param);
  bool clean = true;
  //~ cout << "dsvp call" << endl;
  bkzObj.dSvpReduction(start, param.blockSize, param, clean);
  //~ bkzObj.dumpGSO("gso.log", "DSVP");
  //~ cout << "clean: " << clean << endl;
  return bkzObj.status;
}

// testing
int svpEnum(IntMatrix& b, IntVect& solCoord, bool dual) {
    // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.getRows();

  // Allocates space for vectors and matrices in constructors
  IntMatrix emptyMat;
  MatGSO<Integer, Float> gso(b, emptyMat, emptyMat, GSO_INT_GRAM);
  Float maxDist;
  long maxDistExpo;
  Integer itmp1;

  gso.updateGSO();
  genZeroVect(solCoord, d);
  
  Evaluator<Float>* evaluator;
  evaluator = new FastEvaluator<Float>(d, gso.getMuMatrix(),
            gso.getRMatrix(), EVALMODE_SV);
  evaluator->solCoord.clear();
  
  if (dual) {
    //~ cout << "enumerating dual" << endl;
    maxDist = gso.getRExp(d-1, d-1, maxDistExpo);
    Enumeration::enumerateDual(gso, maxDist, maxDistExpo, *evaluator, 0, d);
  } else {
    vector<double> pruning;
    const vector<Float> emptySubTree;
    maxDist = gso.getRExp(0, 0, maxDistExpo);
    Enumeration::enumerate(gso, maxDist, maxDistExpo, *evaluator, emptySubTree,
            emptySubTree, 0, d, pruning);
  }
  
  if (!evaluator->solCoord.empty()) {
    for (int i = 0; i < d; i++) {
      itmp1.set_f(evaluator->solCoord[i]);
      solCoord[i].add(solCoord[i], itmp1);
    }
  } else {
      solCoord[d-1] = 1;
  }
  
  delete evaluator;
  
  return 0;
}

int dSVPReduction(IntMatrix& b) {
  IntVect x;
  double start = cputime();
  svpEnum(b, x, true);
  cout << "t_enum: " << (cputime() - start)* 0.001 << endl;
  start = cputime();
  
  int d = b.getRows();
  
  // don't want to deal with negativ coefficients
  for (int i = 0; i < d; i++) {
    if (x[i] < 0) {
      x[i].neg(x[i]);
      for (int j = 0; j < b.getCols(); j++) {
        b[i][j].neg(b[i][j]);
      }
    }
  }
  
  int off = 1;
  int k;
	while (off < d) {
		k = d-1;
		while(k - off >= 0) {
			if (x[k] != 0 || x[k-off] != 0) {
        if (x[k] < x[k-off]) {
          x[k].swap(x[k-off]);
          b.swapRows(k, k-off);
        }
        
        while (x[k-off] != 0) {
          while (x[k-off] <= x[k]) {
            x[k].sub(x[k], x[k-off]);
            b[k].sub(b[k-off]);
          }
          
          x[k].swap(x[k-off]);
          b.swapRows(k, k-off);
        }
      }
			k -= 2*off;
    }
		off *= 2;
  }
  
  cout << "t_insert: " << (cputime() - start)* 0.001 << endl;
  return 0;
}
//

FPLLL_END_NAMESPACE
