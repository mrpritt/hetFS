#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <limits>
//#define HOMOGENEOUS
#define RUNNING "running het-FS CrNEH2 ScatterSearch SmallNeigh"
using std::vector;
using std::pair;
void prntV(vector < unsigned short int >&v) {
  for (vector < unsigned short int >::iterator it = v.begin(), _it = v.end(); it != _it; ++it)
    std::cerr << " " << *it;
  std::cerr << std::endl;
}
inline unsigned long int elapsed(bool reset = false) {
  static clock_t start = clock();
  if (reset)
    start = clock();
  return 100.0 * (double) (clock() - start) / (double) CLOCKS_PER_SEC;
}
#define MAXRUNTIME 60000ul
#define MAXRUNITRS 1000000ul
#define MIN1 10
#define MIN2 20
#define BR1 6
#define BR2 4
#define PSZ 60
typedef unsigned int t_time;
typedef unsigned short int t_operation;
typedef unsigned short int t_worker;
typedef unsigned short int t_mach;
const t_time max_time = std::numeric_limits < t_time >::max();
const t_operation no_op = std::numeric_limits < t_operation >::max();
unsigned short int nJobs;
unsigned short int mMach;
#define wWrkr mMach
std::ofstream fout;
vector < t_time > p;
vector < t_mach > op_mach;
vector < t_operation > mach_frst_op;
vector < vector < t_mach > >mach_ops;
typedef struct _solution {
  t_time makespan;
  vector < t_worker > wrkrInMch;
  vector < t_operation > operations;
  vector < size_t > op_pos;
  vector < t_operation > criticalPath;
  _solution():makespan(0), wrkrInMch(mMach), operations(mMach * nJobs + 1), op_pos(mMach * nJobs + 1), criticalPath(2 * mMach + nJobs) {
  };
  t_time Evaluate();
  void AssignRandomWorkers();
  t_time CreateNEH();
  t_time CreateRandom();
  void print();
  void LocalSearch();
} t_solution;
bool compareSol(t_solution S1, t_solution S2) {
  return S1.makespan > S2.makespan;
}
void t_solution::print() {
  fout << "\t{";
  for (auto wim:wrkrInMch)
    fout << " " << wim << ",";
  fout << "}" << std::endl;
  fout << "\t{";
  for (auto op:operations)
    fout << " " << op << ",";
  fout << "}" << std::endl;
  fout << " CritPath";
  for (vector < t_operation >::iterator cp = criticalPath.begin(), _cp = criticalPath.end(); cp != _cp; ++cp) {
    fout << " " << *cp;
  }
  fout << std::endl << elapsed() << " " << makespan << std::endl;
}
void load_operations(std::ifstream & fin) {
  fin >> nJobs >> mMach >> wWrkr;
  op_mach.resize(mMach * nJobs + 1);
  p.resize(wWrkr * nJobs * mMach + wWrkr);
  mach_ops.resize(mMach);
  for (unsigned short int i = 0; i < mMach; ++i) {
    mach_ops[i].reserve(nJobs);
    mach_ops[i].resize(0);
  }
  for (unsigned short int j = 0; j < nJobs; ++j) {
    for (unsigned short int i = 1; i <= mMach; ++i) {
      unsigned short int o = j * mMach + i;
      fin >> op_mach[o];
      mach_ops[i - 1].push_back(o);
      for (unsigned short int k = 0; k < wWrkr; ++k) {
	std::string a;
	fin >> a;
	p[k + o * wWrkr] = (a == "inf") ? max_time : (t_time) stoi(a);
      }
    }
  }
  mach_frst_op.resize(mMach);
  for (unsigned short int i = 1; i <= mMach; ++i)
    mach_frst_op[op_mach[i]] = i;
}
t_time _solution::Evaluate() {
  vector < t_operation > last_op(mMach, 0);
  vector < t_operation > predc;
  vector < t_time > op_compl;
  predc.resize(mMach * nJobs + 1);
  op_compl.resize(mMach * nJobs + 1);
  for (auto op = operations.begin() + 1, _op = operations.end(); op < _op; ++op) {
    unsigned short int mach = op_mach[*op];
    unsigned short int wrkr = wrkrInMch[mach];
    unsigned short int lopm = last_op[mach];
    predc[*op] = (*op % mMach == 1) ? lopm : (op_compl[lopm] < op_compl[*op - 1]) ? *op - 1 : lopm;
    op_compl[*op] = op_compl[predc[*op]] + p[wrkr + *op * wWrkr];
    last_op[mach] = *op;
  }
  makespan = op_compl[last_op.front()];
  t_operation lop = last_op.front();
  for (auto op = last_op.begin() + 1, _op = last_op.end(); op < _op; ++op)
    if (makespan < op_compl[*op]) {
      makespan = op_compl[*op];
      lop = *op;
    }
  criticalPath.resize(0);
  for (; lop > 0; lop = predc[lop])
    criticalPath.push_back(lop);
  return 0;
}
void _solution::AssignRandomWorkers() {
  for (auto it = wrkrInMch.begin(), _it = wrkrInMch.end(); it != _it; ++it)
    *it = it - wrkrInMch.begin();
#ifndef HOMOGENEOUS
  bool invalid = false;
  do {
    invalid = false;
    random_shuffle(wrkrInMch.begin(), wrkrInMch.end());
    for (unsigned short int i = 1; i <= mMach; ++i) {
      if (p[wrkrInMch[op_mach[i]] + i * wWrkr] == max_time) {
	invalid = true;
	break;
      }
    }
  } while (invalid);
#endif
}
t_time _solution::CreateNEH() {
  vector < pair < t_time, t_operation > >totals;
  vector < t_operation > jobs2put, jobs;
  const unsigned mM = mMach;
  const unsigned nJ = nJobs;
  t_time heads[nJ + 1][mM + 1], tails[nJ + 2][mM + 2], rh[nJ + 2][mM + 1];
  totals.resize(nJ);
  jobs2put.resize(nJ);
  jobs.reserve(nJ);
  for (unsigned short int j = 0; j <= nJ; ++j)
    heads[j][0] = tails[j][mM + 1] = rh[j][0];
  tails[nJ + 1][mM + 1] = rh[nJ + 1][0] = 0;
  for (unsigned short int i = 0; i <= mM; ++i)
    heads[0][i] = rh[0][i] = 0;
  AssignRandomWorkers();
  for (unsigned short int j = 0; j < nJ; ++j)
    jobs2put[j] = j;
  random_shuffle(jobs2put.begin(), jobs2put.end());
  jobs.push_back(jobs2put.back());
  for (unsigned short int k = 1; k < nJ; ++k) {
    t_operation job = jobs2put[nJ - k - 1];
    for (unsigned short int i = 1; i <= mM; ++i)
      heads[1][i] = heads[1][i - 1] + p[wrkrInMch[i - 1] + mM * (i + mM * jobs[0])];
    for (unsigned short int j = 2; j <= k; ++j)
      for (unsigned short int i = 1; i <= mM; ++i)
	heads[j][i] = std::max(heads[j][i - 1], heads[j - 1][i])
	  + p[wrkrInMch[i - 1] + mM * (i + mM * jobs[j - 1])];
    for (unsigned short int i = mM; 1 <= i; --i)
      tails[k][i] = tails[k][i + 1] + p[wrkrInMch[i - 1] + mM * (i + mM * jobs[k - 1])];
    for (unsigned short int j = k - 1; 1 <= j; --j)
      for (unsigned short int i = mM; 1 <= i; --i)
	tails[j][i] = std::max(tails[j][i + 1], tails[j + 1][i])
	  + p[wrkrInMch[i - 1] + mM * (i + mM * jobs[j - 1])];
    for (unsigned short int i = 1; i <= mM; ++i)
      rh[1][i] = rh[1][i - 1] + p[wrkrInMch[i - 1] + mM * (i + mM * job)];
    for (unsigned short int j = 2; j <= k + 1; ++j)
      for (unsigned short int i = 1; i <= mM; ++i)
	rh[j][i] = std::max(rh[j][i - 1], heads[j - 1][i])
	  + p[wrkrInMch[i - 1] + mM * (i + mM * job)];
    unsigned short int pos = k + 1;
    t_time Cp = rh[k + 1][mM];
    for (unsigned short int j = 1; j <= k; ++j) {
      t_time Ck = 0;
      for (unsigned short int i = 1; i <= mM; ++i)
	Ck = std::max(Ck, rh[j][i] + tails[j][i]);
      if (Ck < Cp) {
	Cp = Ck;
	pos = j;
      }
    }
    jobs.insert(jobs.begin() + pos - 1, job);
  }
  operations.resize(0);
  operations.push_back(0);
  for (auto jb = jobs.begin(), _jb = jobs.end(); jb < _jb; ++jb) {
    for (t_operation op = *jb * mM + 1, _op = *jb * mM + mM; op <= _op; ++op) {
      op_pos[op] = operations.size();
      operations.push_back(op);
    }
  }
  return Evaluate();
}
t_time _solution::CreateRandom() {
  AssignRandomWorkers();
  vector < t_operation > toSchd;
  toSchd.reserve(2 * nJobs);
  for (unsigned short int i = 1, _i = mMach * nJobs; i <= _i; i += mMach)
    toSchd.push_back(i);
  operations.resize(1);
  const unsigned short int last_op = mMach - 1;
  while (toSchd.size() ) {
    unsigned short int t = rand() % toSchd.size();
    unsigned short int o = toSchd[t];
    op_pos[o] = operations.size();
    operations.push_back(o);
    if (o % mMach)
      {
	toSchd[t]++;
      } else {
      toSchd[t] = toSchd.back();
      toSchd.pop_back();
    }
  }
  return Evaluate();
}
void Position_Preservation_Crossover(vector < t_worker > &P1, vector < t_worker > &P2) {
  vector < unsigned short int >missing;
  vector < unsigned short int >pos, posbk;
  vector < unsigned short int >::iterator it;
  vector < t_worker > bk(P1);
  missing.reserve(wWrkr);
  posbk.reserve(wWrkr);
  pos.reserve(wWrkr);
  posbk.resize(0);
  missing.resize(0);
  missing.resize(P1.size());
  for (unsigned short int i = 0, _i = P1.size(); i < _i; ++i)
    if (P1[i] != P2[i]) {
      missing[i] = 1;
      posbk.push_back(i);
    }
  do {
    pos = posbk;
    for (vector < unsigned short int >::iterator tr = pos.begin(); tr != pos.end(); ++tr)
      missing[P1[*tr] = bk[*tr]] = 1;
    random_shuffle(pos.begin(), pos.end());
    for (vector < unsigned short int >::iterator tr = pos.begin(); tr != pos.end();) {
      unsigned short int t = (rand() % 2) ? P1[*tr] : P2[*tr];
      if (missing[t]) {
	P1[*tr] = t;
	missing[t] = 0;
	*tr = pos.back();
	pos.pop_back();
      } else {
	++tr;
      }
    }
    it = pos.begin();
    for (vector < unsigned short int >::iterator tr = posbk.begin(); tr != posbk.end(); ++tr)
      if (missing[P2[*tr]]) {
	unsigned short int t = P2[*tr];
	if (p[t + mach_frst_op[*it] * wWrkr] == max_time)
	  break;
	P1[*it] = t;
	++it;
      }
  }
  while (it != pos.end());
}
void switch_operations(t_solution & st, vector < t_operation >::iterator ins, vector < t_operation >::iterator end) {
  static vector < t_operation > queue(4 * nJobs);
  queue.resize(0);
  queue.push_back(*ins);
  auto evl = ins + 1;
  for (; evl < end; ++evl) {
    for (auto it = queue.begin(); it < queue.end(); ++it) {
      if ((*evl - 1) / mMach == (*it - 1) / mMach || (*evl) % mMach == (*it) % mMach) {
	queue.push_back(*evl);
	break;
      }
    }
    if (queue.back() != *evl) {
      *ins = *evl;
      st.op_pos[*ins] = ins - st.operations.begin();
      ++ins;
    }
  }
  *ins = *end;
  st.op_pos[*ins] = ins - st.operations.begin();
  ++ins;
  for (auto it = queue.begin(); it < queue.end(); ++it, ++ins) {
    *ins = *it;
    st.op_pos[*ins] = ins - st.operations.begin();
  }
}
void PathRelinking(t_solution SC, t_solution & SG, t_solution & BestSol) {
  t_solution st, ns;
#ifndef HOMOGENEOUS
  Position_Preservation_Crossover(SC.wrkrInMch, SG.wrkrInMch);
#endif
  SC.Evaluate();
  BestSol = SC;
  BestSol.makespan = max_time - 1;
  for (; ;) {
    ns.makespan = max_time;
    for (unsigned short int l = 2, _l = SC.criticalPath.size() - 1; l < _l; ++l) {
      if (SC.criticalPath[l - 1] % mMach == SC.criticalPath[l] % mMach
	  && (SC.criticalPath[l] % mMach != SC.criticalPath[l - 2] % mMach || SC.criticalPath[l] % mMach != SC.criticalPath[l + 1] % mMach)
	  ) {
	if (SG.op_pos[SC.criticalPath[l - 1]] < SG.op_pos[SC.criticalPath[l]]) {
	  st.wrkrInMch = SC.wrkrInMch;
	  st.operations = SC.operations;
	  st.op_pos = SC.op_pos;
	  unsigned short int i = (SC.criticalPath[l] - 1) % mMach;
	  switch_operations(st, st.operations.begin() + st.op_pos[SC.criticalPath[l]], st.operations.begin() + st.op_pos[SC.criticalPath[l - 1]]);
	  st.Evaluate();
	  if (st.makespan < ns.makespan)
	    ns = st;
	}
      }
    }
    if (ns.makespan == max_time) {
      int r = rand() % (SC.criticalPath.size() - 2) + 1;
      for (unsigned short int l = r, _l = SC.criticalPath.size(); l < _l; ++l) {
	if (SC.criticalPath[l - 1] % mMach == SC.criticalPath[l] % mMach) {
	  if (SG.op_pos[SC.criticalPath[l - 1]] < SG.op_pos[SC.criticalPath[l]]) {
	    ns = SC;
	    unsigned short int i = (SC.criticalPath[l] - 1) % mMach;
	    switch_operations(ns, ns.operations.begin() + ns.op_pos[SC.criticalPath[l]], ns.operations.begin() + ns.op_pos[SC.criticalPath[l - 1]]);
	    ns.Evaluate();
	    break;
	  }
	}
      }
      if (ns.makespan == max_time) {
	for (unsigned short int l = 1, _l = r; l < _l; ++l) {
	  if (SC.criticalPath[l - 1] % mMach == SC.criticalPath[l] % mMach) {
	    if (SG.op_pos[SC.criticalPath[l - 1]] < SG.op_pos[SC.criticalPath[l]]) {
	      ns = SC;
	      unsigned short int i = (SC.criticalPath[l] - 1) % mMach;
	      switch_operations(ns, ns.operations.begin() + ns.op_pos[SC.criticalPath[l]], ns.operations.begin() + ns.op_pos[SC.criticalPath[l - 1]]);
	      ns.Evaluate();
	      break;
	    }
	  }
	}
      }
    }
    if (ns.makespan == max_time)
      break;
    SC = ns;
    if (BestSol.makespan > ns.makespan)
      BestSol = ns;
  }
  if (BestSol.makespan == max_time - 1)
    BestSol = SC;
}
void t_solution::LocalSearch() {
  t_solution st, ns;
  vector < t_operation > oppairs;
  for (; ;) {
    oppairs.resize(0);
    oppairs.push_back(no_op);
    for (vector < t_operation >::iterator
	   lm = criticalPath.begin(),
	   cm = criticalPath.begin() + 1, nm = criticalPath.begin() + 2, _m = criticalPath.end(); nm != _m; ++lm, ++cm, ++nm) {
      if (*cm % mMach == *lm % mMach && *cm % mMach != *nm % mMach) {
	if (oppairs.back() != *cm) {
	  oppairs.push_back(*lm);
	  oppairs.push_back(*cm);
	}
      }
      if (*cm % mMach != *lm % mMach && *cm % mMach == *nm % mMach) {
	oppairs.push_back(*cm);
	oppairs.push_back(*nm);
      }
    }
    if (oppairs.size() > 3)
      if (oppairs[1] == criticalPath.front())
	oppairs[2] = oppairs[1] = oppairs[0];
    if (oppairs.back() == criticalPath.back()) {
      oppairs.pop_back();
      oppairs.pop_back();
    }
    ns.makespan = max_time;
    for (; oppairs.back() < no_op;) {
      st.wrkrInMch = wrkrInMch;
      st.operations = operations;
      st.op_pos = op_pos;
      unsigned short int o1 = oppairs.back();
      oppairs.pop_back();
      unsigned short int o2 = oppairs.back();
      oppairs.pop_back();
      switch_operations(st, st.operations.begin() + st.op_pos[o1], st.operations.begin() + st.op_pos[o2]);
      st.Evaluate();
      if (st.makespan < ns.makespan)
	ns = st;
    }
    if (makespan > ns.makespan)
      *this = ns;
    else
      break;
  }
}
unsigned int Distance(const t_solution & SC, const t_solution & SG) {
  unsigned int rslt = 0;
  for (unsigned short int i = 0, _i = mMach; i < _i; ++i) {
    for (auto om = mach_ops[i].begin(), _om = mach_ops[i].end(); om != _om; ++om) {
      for (auto om1 = om + 1; om1 != _om; ++om1) {
	if (SG.op_pos[*om] > SG.op_pos[*om1] && SC.op_pos[*om] < SC.op_pos[*om1])
	  ++rslt;
	else if (SG.op_pos[*om] < SG.op_pos[*om1] && SC.op_pos[*om] > SC.op_pos[*om1])
	  ++rslt;
      }
    }
  }
  for (vector < unsigned short int >::const_iterator wc = SC.wrkrInMch.begin(), _wc = SC.wrkrInMch.end(), wg = SG.wrkrInMch.begin(); wc != _wc;
       ++wc, ++wg)
    if (*wc != *wg)
      ++rslt;
  return rslt;
}
void ScatterSearch() {
  static unsigned int iteration = 0;
  static unsigned int Rsizes = 0;
  static unsigned int randoms = 0;
  static unsigned int btime = 0;
  static unsigned int eachtime = 0;
  iteration = Rsizes = randoms = btime = 0;
  t_solution best;
  t_solution bestNEH;
  vector < t_solution > P, R;
  P.resize(PSZ);
  fout << " initial sols " << RUNNING << std::endl;
  for (vector < t_solution >::iterator iS = P.begin(), _iS = P.end(); iS != _iS; ++iS) {
    (*iS).CreateNEH();
    ++randoms;
  }
  unsigned short int tlmt = 1;
  best = bestNEH = P.back();
  while (P.size()) {
    sort(P.begin(), P.end(), compareSol);
    if (P.back().makespan < best.makespan) {
      best = P.back();
      btime = elapsed();
      fout << " impr " << iteration << " " << randoms << " " << double ((iteration) ? double (Rsizes) /
									double (iteration) : double (0)) << " " << btime << " " << best.
	makespan << std::endl;
      if (tlmt < 4)
	if (tlmt * mMach * nJobs < btime) {
	  fout << " tmpr" << tlmt++ << " " << iteration << " " << randoms << " " << double ((iteration) ? double (Rsizes) /
											    double (iteration) : double (0)) << " " << btime << " " <<
	    best.makespan << std::endl;
	}
      best.print();
    }
    if (elapsed() > eachtime) {
      eachtime += 6000;
      fout << " chk  " << iteration << " " << randoms << " " << double ((iteration) ? double (Rsizes) /
									double (iteration) : double (0)) << " " << elapsed() << " " << btime << " " <<
	best.makespan << std::endl;
    }
    R.resize(0);
    R.push_back(P.back());
    P.pop_back();
    unsigned short int Rsize = BR1;
    {
      unsigned short int MIN = MIN1;
      for (vector < t_solution >::reverse_iterator iS = P.rbegin(), _iS = P.rend(); iS != _iS; ++iS) {
	unsigned short int D = 0;
	for (vector < t_solution >::iterator iR = R.begin(), _iR = R.end(); iR != _iR; ++iR) {
	  D = Distance(*iS, *iR);
	  if (D < MIN)
	    break;
	}
	if (D >= MIN) {
	  R.push_back(*iS);
	  if (R.size() >= Rsize) {
	    if (Rsize == BR1) {
	      Rsize = BR1 + BR2;
	      MIN = MIN2;
	    } else
	      break;
	  }
	}
      }
    }
    Rsizes += (R.size() <= BR1) ? R.size() : BR1;
    if (R.size() <= BR1)
      ++randoms;
    if (Rsize == BR1)
      Rsize = R.size() + BR2;
    P.resize(0);
    while (R.size() < Rsize) {
      R.emplace_back();
      R.back().CreateNEH();
      ++randoms;
      R.back().LocalSearch();
    }
    t_solution t;
    for (vector < t_solution >::iterator iR1 = R.begin(), _iR1 = R.begin() + Rsize - BR2; iR1 != _iR1; ++iR1) {
      for (vector < t_solution >::iterator iR2 = _iR1, _iR2 = R.end(); iR2 != _iR2; ++iR2) {
	P.emplace_back();
	PathRelinking(*iR2, *iR1, P.back());
	P.back().LocalSearch();
      }
    }
    if (++iteration > MAXRUNITRS)
      break;
    if (elapsed() > MAXRUNTIME)
      break;
  }
  sort(P.begin(), P.end(), compareSol);
  if (P.back().makespan < best.makespan) {
    best = P.back();
    btime = elapsed();
  }
  fout << " =*= BEST NEH =*= " << std::endl;
  bestNEH.print();
  fout << " =*= BEST =*= " << std::endl;
  best.print();
  fout << "Rslt " << iteration << " " << btime << " " << elapsed() << " " << best.makespan << " " << bestNEH.makespan << std::endl;
}
void RandomSearch() {
  static unsigned int iteration = 0;
  static unsigned int Rsizes = 0;
  static unsigned int randoms = 0;
  static unsigned int btime = 0;
  static unsigned int bNtime = 0;
  iteration = Rsizes = randoms = btime = 0;
  t_solution best, bestNEH, bnow;
  vector < t_solution > P, R;
  P.resize(PSZ);
  best.CreateNEH();
  bestNEH = best;
  best.LocalSearch();
  for (;;) {
    bnow.CreateNEH();
    if (bnow.makespan < bestNEH.makespan) {
      bestNEH = bnow;
      bNtime = elapsed();
      fout << " impN " << iteration << " " << randoms << " " << 0 << " " << bNtime << " " << best.makespan << std::endl;
      if (iteration > 1000)
	best.print();
    }
    bnow.LocalSearch();
    if (bnow.makespan < best.makespan) {
      best = bnow;
      btime = elapsed();
      fout << " impr " << iteration << " " << randoms << " " << 0 << " " << btime << " " << best.makespan << std::endl;
      if (iteration > 1000)
	best.print();
    }
    ++iteration;
    if (elapsed() > MAXRUNTIME)
      break;
  }
  fout << " =*= BEST NEH =*= " << std::endl;
  bestNEH.print();
  fout << " =*= BEST =*= " << std::endl;
  best.print();
  fout << "Rslt " << iteration << " " << btime << " " << elapsed() << " " << best.makespan << " " << bestNEH.makespan << std::endl;
}
int main(int argc, char **argv) {
  int rank, size, namelen;
  int rtst = 1;
  unsigned int seed;
  std::vector < std::string > files;
  files.reserve(50);
  srand(time(0));
  size = 1;
  if (argc > 2) {
    rank = std::stoi(std::string(argv[1]));
    size = std::stoi(std::string(argv[2]));
    seed = rand() + rank * 10000 + size * 100;
    rtst = (argc > 3) ? std::stoi(std::string(argv[3])) : 0;
  } else {
    rank = 0;
    seed = rand();
    size = 1;
    rtst = 0;
  }
  std::ifstream fin("jobs.txt");
  std::string file;
  file.reserve(20);
  for (short int i = -1; i < rank; ++i)
    fin >> file;
  for (;;) {
    if (file.size() > 3)
      files.push_back(std::move(file));
    else
      break;
    file.reserve(20);
    for (short int i = 0; i < size; ++i)
      fin >> file;
  }
  fin.close();
  if (file.size() > 3)
    std::cerr << " Finish jobs.txt with 'eof'." << std::endl;
  file.reserve(40);
  for (auto fl:files)
    std::cout << fl << " ";
  std::cout << std::endl;
  for (auto fl = files.begin(); fl < files.end(); ++fl) {
    file.assign("../instances/");
#ifdef HOMOGENEOUS
    file.assign("../instancesH/");
#endif
    file.append(*fl);
    fin.open(file);
    file.assign(*fl);
    load_operations(fin);
    fin.close();
    *(file.end() - 3) = 'o';
    srand(seed);
    std::cout << "initseed " << seed << std::endl;
    for (short int i = rtst + 1; i <= rtst + 10 ; ++i) {
      *(file.end() - 2) = '0' + (i / 10);
      *(file.end() - 1) = '0' + (i % 10);
      fout.open(file);
      seed = rand();
      fout << "seed " << seed << " " << rank << " " << size << std::endl;
      srand(seed);
      elapsed( true);
      ScatterSearch();
      fout.close();
    }
  }
  return 0;
}
void test() {
  t_solution S;
  vector < vector < t_operation > >ttt;
  ttt.resize(mMach);
  S.CreateRandom();
  S.print();
  S.wrkrInMch = {
    0, 0, 1, 4, 3, 7, 6, 5, 2};
  ttt[0] = {
    0, 9, 1, 49, 17, 25, 57, 33, 41,};
  ttt[1] = {
    0, 10, 2, 50, 18, 26, 58, 34, 42,};
  ttt[2] = {
    0, 11, 3, 51, 19, 27, 59, 35, 43,};
  ttt[3] = {
    0, 12, 4, 52, 20, 28, 60, 36, 44,};
  ttt[4] = {
    0, 13, 5, 53, 21, 29, 61, 37, 45,};
  ttt[5] = {
    0, 14, 6, 54, 22, 30, 62, 46, 38,};
  ttt[6] = {
    0, 15, 7, 55, 23, 31, 63, 47, 39,};
  ttt[7] = {
    0, 16, 8, 56, 24, 64, 32, 48, 40,};
  for (int i = 0; i < mMach; ++i) {
    for (int j = 1; j <= nJobs; ++j)
      S.operations[i * nJobs + j] = ttt[i][j];
  }
  S.Evaluate();
  S.print();
}
