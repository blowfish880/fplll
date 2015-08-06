/* Copyright (C) 2011 Xavier Pujol
   (C) 2014 Martin Albrecht.

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

#include <iomanip>

/* Template source file */
#include "bkz.h"
#include "enumerate.h"
#include <iomanip>

FPLLL_BEGIN_NAMESPACE

template<class FT>
bool BKZAutoAbort<FT>::testAbort(double scale, int maxNoDec) {
  double newSlope = -getCurrentSlope(m, startRow, numRows);
  if (noDec == -1 || newSlope < scale*oldSlope)
    noDec = 0;
  else
    noDec++;
  oldSlope = min(oldSlope, newSlope);
  return noDec >= maxNoDec;
}

template<class FT>
BKZReduction<FT>::BKZReduction(MatGSO<Integer, FT>& m,
        LLLReduction<Integer, FT>& lllObj, const BKZParam& param) :
  status(RED_SUCCESS), param(param), m(m), lllObj(lllObj)
{
  for (numRows = m.d; numRows > 0 && m.b[numRows - 1].is_zero(); numRows--) {}
  this->delta = param.delta;
}

template<class FT>
BKZReduction<FT>::~BKZReduction() {
}

template<class FT>
static double getCurrentSlope(MatGSO<Integer, FT>& m, int startRow, int stopRow) {
  FT f, logF;
  long expo;
  vector<double> x;
  x.resize(stopRow);
  for (int i = startRow; i < stopRow; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    x[i] = logF.get_d() + expo * log(2.0);
  }
  int n = stopRow - startRow;
  double iMean = (n - 1) * 0.5 + startRow, xMean = 0, v1 = 0, v2 = 0;
  for (int i = startRow; i < stopRow; i++) {
    xMean += x[i];
  }
  xMean /= n;
  for (int i = startRow; i < stopRow; i++) {
    v1 += (i - iMean) * (x[i] - xMean);
    v2 += (i - iMean) * (i - iMean);
  }
  return v1 / v2;
}

template<class FT>
void BKZReduction<FT>::getSubDet(double& subdet, int start, int end) {
  subdet = 0.0;
  FT f, logF;
  long expo;
  for(int i=start; i<end; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    
    subdet += logF.get_d() + expo * log(2.0);
  }
}

template<class FT>
void BKZReduction<FT>::setRollback() {
  rb_point = IntMatrix(m.b);
}

template<class FT>
void BKZReduction<FT>::rollback() {
  m.rowOpBegin(0, numRows);
  m.b = rb_point;
  m.rowOpEnd(0, numRows);
  m.updateGSO();
}

template<class FT>
void BKZReduction<FT>::updateSldPotential() {
  sldPot = 0.0;
  int p = numRows/param.blockSize;
  FT f, logF;
  long expo;
  for(int i=0; i<(p*param.blockSize); i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    
    sldPot += (p-(i/param.blockSize))*(logF.get_d() + expo * log(2.0));
  }
}

template<class FT>
bool BKZReduction<FT>::dSvpReduction(int kappa, int blockSize, const BKZParam &par, bool& clean) {
  long maxDistExpo;
  int lllEnd = (par.flags & BKZ_BOUNDED_LLL) ? kappa + blockSize : m.d;
  //~ cout << "LLL call from " << kappa << " to " << kappa + blockSize << endl;
  if (!lllObj.lll(kappa, kappa, lllEnd)) {
    return setStatus(lllObj.status);
  }
  if (lllObj.nSwaps > 0) {
    clean = false;
  }
  
  //~ cout << "LLL done. preparing enum" << endl;
  
  if (!localPP(par, kappa, blockSize, clean, true)) {
    return false;
  }
  
  m.updateGSO();
  maxDist = m.getRExp(kappa + blockSize - 1, kappa + blockSize - 1, maxDistExpo);
  deltaMaxDist.div(maxDist, delta);
  //~ cout << "maxDist " <<  maxDist << endl;
  //~ cout << "deltaMaxDist " <<  deltaMaxDist << endl;
  vector<FT>& x = evaluator.solCoord;
  x.clear();
  //~ cout << "enum call" << endl;
  //~ Enumeration::enumerateDual(m, maxDist, maxDistExpo, evaluator, kappa, kappa + blockSize, par.pruning);
  Enumeration::enumerate(m, maxDist, maxDistExpo, evaluator, emptySubTree,
            emptySubTree, kappa, kappa + blockSize, par.pruning, true);
  
  if (maxDist <= deltaMaxDist) {
    // basis already delta-dSVP reduced
    //~ cout << maxDist << " <= " << deltaMaxDist << endl;
    return true;
  }
  
  //~ cout << "non trivial solution found. inserting: " << x << endl;
  // we don't have to worry about the vector already being part of the basis
  // the gcd procedure will take care of it automatically and very efficiently
  int d = blockSize;
  m.rowOpBegin(kappa, kappa + blockSize);
  // don't want to deal with negativ coefficients
  for (int i = 0; i < d; i++) {
    if (x[i] < 0) {
      x[i].neg(x[i]);
      for (int j = 0; j < m.b.getCols(); j++) {
        m.b[i + kappa][j].neg(m.b[i + kappa][j]);
      }
    }
  }
  //~ cout << "enum done. inserting..." << endl;
  // tree based gcd computation on x, performing operations also on b
  int off = 1;
  int k;
	while (off < d) {
		k = d - 1;
		while(k - off >= 0) {
			if (!(x[k].is_zero() && x[k - off].is_zero())) {
        if (x[k] < x[k - off]) {
          x[k].swap(x[k - off]);
          m.b.swapRows(kappa + k, kappa + k - off);
        }
        
        while (!x[k - off].is_zero()) {
          while (x[k - off] <= x[k]) {
            x[k].sub(x[k], x[k - off]);
            m.b[kappa + k].sub(m.b[kappa + k - off]);
          }
          
          x[k].swap(x[k - off]);
          m.b.swapRows(kappa + k, kappa + k - off);
        }
      }
			k -= 2 * off;
    }
		off *= 2;
  }
  //~ cout << "insertion loop done" << endl;
  m.rowOpEnd(kappa, kappa + blockSize);
  clean = false;
  if (!lllObj.lll(kappa, kappa, kappa + blockSize)) {
    return setStatus(lllObj.status);
  }
  return true;
}

template<class FT>
bool BKZReduction<FT>::localPP(const BKZParam &par, int kappa, int blockSize, bool& clean, bool dual) {
  ppCputimeStart = cputime();
  const BKZParam *preproc = par.preprocessing;
  if (preproc && preproc->blockSize < blockSize && preproc->blockSize > 2) {
    int dummyKappaMax = numRows;
    BKZAutoAbort<FT> autoAbort(m, kappa + blockSize, kappa);
    double cputimeStart2 = cputime();

    for(int i=0; ; i++) {
      if ((preproc->flags & BKZ_MAX_LOOPS) && i >= preproc->maxLoops) break;
      if ((preproc->flags & BKZ_MAX_TIME) && (cputime() - cputimeStart2) * 0.001 >= preproc->maxTime) break;
      if (autoAbort.testAbort(preproc->autoAbort_scale, preproc->autoAbort_maxNoDec)) break;

      bool clean2 = true;
      if (dual) {
        if (!dbkzLoop(i, dummyKappaMax, *preproc, kappa, kappa + blockSize, clean2))
          return false;
      } else {
        if (!bkzLoop(i, dummyKappaMax, *preproc, kappa, kappa + blockSize, clean2))
          return false;
      }

      if(clean2)
        break;
      else
        clean = clean2;
    }
  }
  ppCputime += (cputime() - ppCputimeStart);
  
  return true;
}

template<class FT>
bool BKZReduction<FT>::svpReduction(int kappa, int blockSize, const BKZParam &par, bool& clean) {
  long maxDistExpo;

  int lllStart = (par.flags & BKZ_BOUNDED_LLL) ? kappa : 0;

  if (!lllObj.lll(lllStart, kappa, kappa + blockSize)) {
    return setStatus(lllObj.status);
  }
  if (lllObj.nSwaps > 0) {
    clean = false;
  }
  
  if (!localPP(par, kappa, blockSize, clean, false)) {
    return false;
  }
  
  m.updateGSO();
  maxDist = m.getRExp(kappa, kappa, maxDistExpo);
  deltaMaxDist.mul(delta, maxDist);
  vector<FT>& solCoord = evaluator.solCoord;
  solCoord.clear();
  Enumeration::enumerate(m, maxDist, maxDistExpo, evaluator, emptySubTree,
            emptySubTree, kappa, kappa + blockSize, par.pruning);
  if (solCoord.empty()) {
    return setStatus(RED_ENUM_FAILURE);
  }

  // Is it already in the basis ?
  int nzVectors = 0, iVector = -1;
  for (int i = 0; i < blockSize; i++) {
    if (!solCoord[i].is_zero()) {
      nzVectors++;
      if (iVector == -1 && (solCoord[i].get_d() == 1 || solCoord[i].get_d() == -1))
        iVector = i;
    }
  }
  FPLLL_DEBUG_CHECK(nzVectors > 0);

  if (maxDist >= deltaMaxDist) {
    return true; // Do nothing
  }

  if (nzVectors == 1) {
    // Yes, it is another vector
    FPLLL_DEBUG_CHECK(iVector != -1 && iVector != 0);
    m.moveRow(kappa + iVector, kappa);
    if (!lllObj.sizeReduction(kappa, kappa + 1))
      return setStatus(lllObj.status);
  }
  else {
    // No, general case
    int d = m.d;
    m.createRow();
    m.rowOpBegin(d, d + 1);
    for (int i = 0; i < blockSize; i++) {
      m.row_addmul(d, kappa + i, solCoord[i]);
    }
    m.rowOpEnd(d, d + 1);
    m.moveRow(d, kappa);
    if (!lllObj.lll(kappa, kappa, kappa + blockSize + 1))
      return setStatus(lllObj.status);
    FPLLL_DEBUG_CHECK(m.b[kappa + blockSize].is_zero());
    m.moveRow(kappa + blockSize, d);
    m.removeLastRow();
  }
  clean = false;
  return true;
}

template<class FT>
bool BKZReduction<FT>::bkzLoop(const int loop, int& kappaMax, const BKZParam &par, int minRow, int maxRow, bool& clean) {
  for (int kappa = minRow; kappa < maxRow-1; kappa++) {
    // SVP-reduces a block
    int blockSize = min(par.blockSize, maxRow - kappa);
    if (!svpReduction(kappa, blockSize, par, clean)) return false;
    if ((par.flags & BKZ_VERBOSE) && kappaMax < kappa && clean) {
      cerr << "Block [1-" << setw(4) << kappa + 1 << "] BKZ-" << setw(0) << par.blockSize << " reduced for the first time" << endl;
      kappaMax = kappa;
    }
  }

  if (par.flags & BKZ_VERBOSE) {
    FT r0;
    Float fr0;
    long expo;
    r0 = m.getRExp(minRow, minRow, expo);
    fr0 = r0.get_d();
    fr0.mul_2si(fr0, expo);
    cerr << "End of BKZ loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s";
    cerr << ", r_" << minRow << " = " << fr0;
    cerr << ", slope = " << std::setw( 9 ) << std::setprecision( 6 ) << getCurrentSlope(m, minRow, maxRow) << endl;
  }
  if (par.flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "End of BKZ loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s, ";
    prefix << std::setw(20) << Enumeration::nodes << ")";
    dumpGSO(par.dumpGSOFilename, prefix.str());
  }

  return true;
}

template<class FT>
bool BKZReduction<FT>::dbkzLoop(const int loop, int& kappaMax, const BKZParam &par, int minRow, int maxRow, bool& clean) {
  m.updateGSO();
  for (int kappa = maxRow-par.blockSize; kappa >= minRow; kappa--) {
    // dSVP-reduces a block
    if (!dSvpReduction(kappa, par.blockSize, par, clean)) return false;
  }
  
  for (int blockSize = par.blockSize - 1; blockSize > 1; blockSize--) {
    if (!dSvpReduction(minRow, blockSize, par, clean)) return false;
  }
  
  if (par.flags & BKZ_VERBOSE) {
    FT r0;
    Float fr0;
    long expo;
    r0 = m.getRExp(minRow, minRow, expo);
    fr0 = r0.get_d();
    fr0.mul_2si(fr0, expo);
    cerr << "End of DBKZ loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s";
    cerr << ", r_" << minRow << " = " << fr0;
    cerr << ", slope = " << std::setw( 9 ) << std::setprecision( 6 ) << getCurrentSlope(m, minRow, maxRow) << endl;
  }
  if (par.flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "End of DBKZ loop " << std::setw(3) << loop;
    prefix << " (" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s,";
    prefix << std::setw(20) << Enumeration::nodes << ")";
    dumpGSO(par.dumpGSOFilename, prefix.str());
  }

  return true;
}

template<class FT>
bool BKZReduction<FT>::rnkLoop(const int loop, const BKZParam &par, int minRow, int maxRow, bool& clean) {
  //~ double subDetLp;
  int dummyKappaMax = numRows;
  FT logDelta; logDelta.log(delta, GMP_RNDU);
  bool clean2;
  int lp;
  
  for (int kappa = minRow; kappa < maxRow - par.blockSize - 1; kappa += par.blockSize) {
    setRollback();
    getSubDet(subDetOld, kappa, kappa + par.blockSize);
    //~ subDet = subDetOld;
    lp = 0;
    //~ do {
      //~ subDetLp = subDet;
      //~ 
      //~ if (!dbkzLoop(lp, dummyKappaMax, par, kappa, min(kappa + 2*par.blockSize, maxRow), clean2))
        //~ return false;
      //~ 
      //~ if (!bkzLoop(lp, dummyKappaMax, par, kappa, min(kappa + 2*par.blockSize, maxRow), clean2))
        //~ return false;
      //~ 
      //~ getSubDet(subDet, kappa, kappa + par.blockSize);
      //~ 
      //~ if (subDet <= subDetOld + logDelta.get_d() && subDet < subDetLp) {
        //~ setRollback();
        //~ clean = false;
      //~ }
      //~ lp++;
    //~ } while (subDet < subDetLp);
    
    int cnt = -1;
    double oldSlope, newSlope;
    while (cnt < 5) {
      if (!dbkzLoop(lp, dummyKappaMax, par, kappa, min(kappa + 2*par.blockSize, maxRow), clean2))
        return false;
      
      if (!bkzLoop(lp, dummyKappaMax, par, kappa, min(kappa + 2*par.blockSize, maxRow), clean2))
        return false;
        
      newSlope = -getCurrentSlope(m, kappa, min(kappa + 2*par.blockSize, maxRow));
      if (cnt < 0 || newSlope < oldSlope) {
        cnt = 0;
      } else {
        cnt++;
      }
      oldSlope = min(oldSlope, newSlope);
      lp++;
    }
    
    getSubDet(subDet, kappa, kappa + par.blockSize);
    if (subDet > subDetOld + logDelta.get_d()) {
      rollback();
    } else {
      clean = false;
    }
  }
  
  if (par.flags & BKZ_VERBOSE) {
    FT r0;
    Float fr0;
    long expo;
    r0 = m.getRExp(minRow, minRow, expo);
    fr0 = r0.get_d();
    fr0.mul_2si(fr0, expo);
    cerr << "End of RNK loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s";
    cerr << ", r_" << minRow << " = " << fr0;
    cerr << ", slope = " << std::setw( 9 ) << std::setprecision( 6 ) << getCurrentSlope(m, minRow, maxRow) << endl;
  }
  
  if (par.flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "End of RNK loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s, ";
    prefix << std::setw(20) << Enumeration::nodes << ")";
    dumpGSO(par.dumpGSOFilename, prefix.str());
  }
  
  return true;
}

template<class FT>
bool BKZReduction<FT>::sldLoop(const int loop, const BKZParam &par, int minRow, int maxRow, bool& clean) {
  bool clean_inner;
  //~ std::ostringstream preLLL, preSVP, preDSVP;
  //~ preLLL << "LLL ";
  //~ preSVP << "SVP ";
  //~ preDSVP << "DSVP ";
  
  do {
    clean_inner = true;
    
    if (!lllObj.lll(0, minRow, maxRow)) {
      return setStatus(lllObj.status);
    }
    if (lllObj.nSwaps > 0) {
      //~ cout << "LLL changed basis" << endl;
      clean_inner = false;
    }
    
    //~ dumpGSO(par.dumpGSOFilename, preLLL.str());
    
    for (int kappa = minRow; kappa < maxRow - 1; kappa += par.blockSize) {
      int blockSize = min(par.blockSize, maxRow - kappa);
      if (!svpReduction(kappa, blockSize, par, clean_inner)) {
        //~ cout << "SVP failed for block " << kappa << ":" << (kappa + blockSize) << endl;
        return false;
      }
      svpCalls++;
      m.updateGSO();
    }
    
    //~ dumpGSO(par.dumpGSOFilename, preSVP.str());
    
    if (!clean_inner) {
      //~ cout << "SVP changed basis" << endl;
      clean = false;
    }
  } while (!clean_inner);
  
  for (int kappa = minRow + 1; kappa < maxRow - par.blockSize; kappa += par.blockSize) {
    if (!dSvpReduction(kappa, par.blockSize, par, clean)) return false;
    svpCalls++;
    m.updateGSO();
  }
  //~ dumpGSO(par.dumpGSOFilename, preDSVP.str());
  
  updateSldPotential();
  clean = (clean || (sldPot >= sldPotOld));
  sldPotOld = sldPot;
  
  if (par.flags & BKZ_VERBOSE) {
    FT r0;
    Float fr0;
    long expo;
    r0 = m.getRExp(minRow, minRow, expo);
    fr0 = r0.get_d();
    fr0.mul_2si(fr0, expo);
    cerr << "End of SLD loop " << std::setw(4) << loop << ", time = " << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s";
    cerr << ", r_" << minRow << " = " << fr0;
    cerr << ", slope = " << std::setw( 9 ) << std::setprecision( 6 ) << getCurrentSlope(m, minRow, maxRow) << endl;
  }
  
  if (par.flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "End of SLD loop " << std::setw(4) << loop;
    prefix << " (" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart) * 0.001 << "s, ";
    prefix << std::setw(9) << svpCalls << ")";
    dumpGSO(par.dumpGSOFilename, prefix.str());
  }
  
  return true;
}

template<class FT>
bool BKZReduction<FT>::bkz() {
  if (param.red_method != RED_BKZ &&
      param.red_method != RED_SLD &&
      param.red_method != RED_DUAL &&
      param.red_method != RED_DBKZ &&
      param.red_method != RED_RNK) {
    FPLLL_ABORT("Unsupported reduction method!");
  }
  
  int flags = param.flags;
  int finalStatus = RED_SUCCESS;

  if (flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "Input";
    dumpGSO(param.dumpGSOFilename, prefix.str(), false);
  }

  if (param.blockSize < 2)
    return setStatus(RED_SUCCESS);

  int kappaMax = 0;
  int iLoop =0;
  BKZAutoAbort<FT> autoAbort(m, numRows);
  if ( (param.red_method == RED_DBKZ) && !(flags & BKZ_AUTO_ABORT) ) {
    cerr << "Warning: need auto abort for dual BKZ. Turning it on!" << endl;
    flags |= BKZ_AUTO_ABORT;
  }
  
  if (flags & BKZ_VERBOSE) {
    if (param.red_method == RED_SLD) {
      cerr << "Entering SLD:" << endl;
    } else if (param.red_method == RED_DBKZ) {
      cerr << "Entering DBKZ:" << endl;
    } else if (param.red_method == RED_DUAL) {
      cerr << "Entering DUAL:" << endl;
    } else if (param.red_method == RED_BKZ){
      cerr << "Entering BKZ:" << endl;
    } else if (param.red_method == RED_RNK){
      cerr << "Entering RNK:" << endl;
    }
    printParams(param, cerr);
    cerr << endl;
  }
  
  if (param.red_method == RED_SLD) {
    updateSldPotential();
    sldPotOld = sldPot;
  }
  
  cputimeStart = cputime();
  ppCputime = 0;
  svpCalls = 0;
  Enumeration::nodes = 0;

  m.discoverAllRows();

  for (iLoop = 0;; iLoop++) {
    if ((flags & BKZ_MAX_LOOPS) && iLoop >= param.maxLoops) {
      finalStatus = RED_BKZ_LOOPS_LIMIT;
      break;
    }
    if ((flags & BKZ_MAX_TIME) && (cputime() - cputimeStart) * 0.001 >= param.maxTime) {
      finalStatus = RED_BKZ_TIME_LIMIT;
      break;
    }
    if ((flags & BKZ_AUTO_ABORT) && autoAbort.testAbort(param.autoAbort_scale, param.autoAbort_maxNoDec)) break;
    bool clean = true;
    if (param.red_method == RED_SLD) {
      if (!sldLoop(iLoop, param, 0, numRows, clean)) return false;
    } else if (param.red_method == RED_RNK) {
      if (!rnkLoop(iLoop, param, 0, numRows, clean)) return false;
    } else {
      if (param.red_method == RED_DBKZ || param.red_method == RED_DUAL) {
        if (!dbkzLoop(iLoop, kappaMax, param, 0, numRows, clean)) return false;
      }
      
      if (param.red_method == RED_DBKZ || param.red_method == RED_BKZ) {
        if (!bkzLoop(iLoop, kappaMax, param, 0, numRows, clean)) return false;
      }
    }
    if (clean || param.blockSize >= numRows) break;
  }
  if (flags & BKZ_DUMP_GSO) {
    std::ostringstream prefix;
    prefix << "Output ";
    prefix << " (" << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (cputime() - cputimeStart)* 0.001 << "s, ";
    prefix << std::fixed << std::setw( 9 ) << std::setprecision( 3 ) << (ppCputime * 0.001) << "s)";
    
    dumpGSO(param.dumpGSOFilename, prefix.str());
  }
  return setStatus(finalStatus);
}

template<class FT>
void BKZReduction<FT>::printParams(const BKZParam &param, ostream &out) {
  out << "blocksize: " << std::setw(3) << param.blockSize << ", ";
  out << "flags: 0x" << std::setw(4) << setfill('0') << std::hex << param.flags << ", " << std::dec << std::setfill(' ');
  out << "maxLoops: " << std::setw(3) << param.maxLoops << ", ";
  out << "maxTime: " << std::setw(0) << std::fixed << std::setprecision( 1 ) << param.maxTime << ", ";
  if (param.flags & BKZ_AUTO_ABORT) {
    out << "autoAbort: (" << std::setw(0) << std::fixed << std::setprecision( 4 ) << param.autoAbort_scale;
    out << ", " << std::setw(2) << param.autoAbort_maxNoDec << "), ";
  } else {
    out << "autoAbort: (     -,  -), ";
  }
  out << "pruning: ";
  if (!param.pruning.empty())
    out << 1;
  else
    out << 0;
  out << endl;

  if (param.preprocessing)
    printParams(*param.preprocessing, out);
}

template<class FT>
bool BKZReduction<FT>::setStatus(int newStatus) {
  status = newStatus;
  if (param.flags & BKZ_VERBOSE) {
    if (status == RED_SUCCESS)
      cerr << "End of BKZ: success" << endl;
    else
      cerr << "End of BKZ: failure: " << RED_STATUS_STR[status] << endl;
  }
  return status == RED_SUCCESS;
}

template<class FT>
void BKZReduction<FT>::dumpGSO(const std::string filename, const std::string prefix, bool append) {
  ofstream dump;
  if (append)
    dump.open(filename.c_str(), std::ios_base::app);
  else
    dump.open(filename.c_str());
  dump << std::setw(4) << prefix << ": ";
  FT f, logF;
  long expo;
  for(int i=0; i<numRows; i++) {
    m.updateGSORow(i);
    f = m.getRExp(i, i, expo);
    logF.log(f, GMP_RNDU);
    
    if (isnan(logF.get_d())) {
      cerr << m.getMuMatrix() << endl;
      cerr << endl;
      cerr << m.getRMatrix() << endl;
      
      cout << i << ": " << f << endl;
      cout << m.updateGSO() << endl;
      cout << i << ": " << f << endl;
      
      dump << std::endl;
      dump.close();
      FPLLL_ABORT("NaN detected!");
    }
    
    dump << std::setprecision(8) << logF.get_d() + expo * log(2.0) << " ";
  }
  dump << std::endl;
  dump.close();
}


FPLLL_END_NAMESPACE
