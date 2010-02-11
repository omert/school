#include <iomanip>
#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include "Matrix.hpp"

using namespace std;

struct GF27 : public map<int, int> {

    void normalize();
    void print() const;
};

GF27
operator + (const GF27& a, const GF27& b)
{
    GF27 ret;
    for (GF27::const_iterator it = a.begin(); it != a.end(); ++it)
	ret[it->first] += it->second;
    for (GF27::const_iterator it = b.begin(); it != b.end(); ++it)
	ret[it->first] += it->second;
    ret.normalize();
    return ret;
}

GF27
operator * (const GF27& a, const GF27& b)
{
    GF27 ret;
    for (GF27::const_iterator ita = a.begin(); ita != a.end(); ++ita)
	for (GF27::const_iterator itb = b.begin(); itb != b.end(); ++itb)
	    ret[ita->first + itb->first] += ita->second * itb->second;
    ret.normalize();
    return ret;
}

void
GF27::print() const
{
    for (const_reverse_iterator it = rbegin(); it != rend(); ++it){
	if (it != rbegin()){
	    if (it->second > 0)
		cout << '+';
	    else
		cout << '-';
	}
	if (it->first > 0){
	    if (it->second > 1)
		cout << it->second;
	    cout << "a";
	    if (it->first > 1)
		cout << "^{" << it->first << "}";
	}
	else
	    cout << it->second;

    }
}

void
GF27::normalize()
{
    iterator itNext = begin();
    for (iterator it = begin(); it != end(); it = itNext){
	++itNext;
	it->second %= 3;
	it->second += 3;
	it->second %= 3;
	if (it->second == 0)
	    erase(it);
    }
    if (size() && rbegin()->first > 2){
	int exp = rbegin()->first;
	int coeff = rbegin()->second;
	erase(exp);
	GF27 a;
	a[exp-3] = coeff;
	GF27 b;
	b[0] = -1;
	b[1] = -2;
	*this = *this + a*b;
	normalize();
    }
}



void
rootTable()
{
    for (size_t i = 0; i < 26; ++i){
	GF27 x;
	x[i] = 1;
	x.normalize();
	cout << "$a^{" << i << "}$" << " & ";
	cout << '$';
	x.print();
	cout << '$';
	cout << "\\\\" << endl;
    }
	
}

void
minpolyTable()
{
    GF27 zero;
    GF27 one;
    one[0] = 1;
    GF27 two;
    two[0] = 2;
    for (size_t i = 0; i < 26; ++i){
	GF27 x;
	x[i] = 1;
	x.normalize();
	cout << "$a^{" << i << "}$" << " & ";
	if (x + one == zero)
	    cout << "$" << "x-2" << "$" << " & " << "(i)";
	if (x + two == zero)
	    cout << "$" << "x-1" << "$" << " & " << "(i)";
	if (x * x * x + two * x + one == zero)
	    cout << "$" << "x^3+2x+1" << "$"  << " & " << "(ii)";
	if (x * x * x + two * x + two == zero)
	    cout << "$" << "x^3+2x+2" << "$" << " & " << "(vi)";
	if (x * x * x + x * x + two == zero)
	    cout << "$" << "x^3+x^2+2" << "$"  << " & " << "(vi)";
	if (x * x * x + x * x + x + two == zero)
	    cout << "$" << "x^3+x^2+x+2" << "$" << " & " << "(v)";
	if (x * x * x + x * x + two * x + one == zero)
	    cout << "$" << "x^3+x^2+2x+1" << "$"  << " & " << "(iv)";
	if (x * x * x + two * x * x + one == zero)
	    cout << "$" << "x^3+2x^2+1" << "$" << " & " << "(iii)";
	if (x * x * x + two * x * x + x + one == zero)
	    cout << "$" << "x^3+2x^2+x+1" << "$"  << " & " << "(iv)";
	if (x * x * x + two * x * x + two * x + two == zero)
	    cout << "$" << "x^3+2x^2+2x+2" << "$" << " & " << "(v)";
	cout << "\\\\" << endl;
	
    }
}


struct Poly : public deque<int> {
    Poly() {}
    Poly(const string& s) {
	for (size_t i = 0; i < s.size(); ++i)
	    push_front(s[i] - '0');
	normalize();
    }

    void normalize();
};


ostream&
operator << (ostream& os, const Poly& p)
{
    for (size_t i = p.size(); i > 0; --i)
	if (p[i - 1] != 0){
	    if (i < p.size())
		os << " + ";
	    if (i == 1)
		os << p[i - 1];
	    else{
		if (p[i - 1] != 1)
		    os << p[i - 1];
		if (i == 2)
		    os << "x";
		else
		    os << "x^" << i - 1;
	    }
	}
    return os;
}

void
Poly::normalize()
{
    for (size_t i = 0; i < size(); ++i){
	(*this)[i] %= 3;
	(*this)[i] += 3;
	(*this)[i] %= 3;
    }
    while (size() && back() == 0)
	pop_back();
}

Poly
operator * (const Poly& p, const Poly& q)
{
    Poly res;
    res.resize(p.size() * q.size());
    for (size_t i = 0; i < p.size(); ++i)
	for (size_t j = 0; j < q.size(); ++j)
	    res[i + j] += p[i] * q[j];
	
    res.normalize();
    return res;
}

void
divide(Poly p1, const Poly& p2, Poly& q, Poly& r)
{
    q.resize(p1.size());
    while (p1.size() >= p2.size()){
	int f = p1.back() * p2.back();
	int e = p1.size() - p2.size();
	q[e] = f;
	for (size_t i = 0; i < p2.size(); ++i)
	    p1[i + e] -= f * p2[i];
	p1.normalize();
    }
    r = p1;
}

void
polyMath(int argc, char* argv[])
{
    if (argc < 3){
	cout << "need two args" << endl;
	return;
    }
    string p1s(argv[1]);
    string p2s(argv[2]);
    Poly p1(p1s);
    Poly p2(p2s);
    Poly q;
    Poly r;
    divide(p1, p2, q, r);
    cout << p1 << " divided by " << p2 << " = " 
	 << q << " with reminder " << r << endl;
    cout << p1 << " times " << p2 << " = " << p1 * p2 << endl;
}

void
cubicPolys()
{
    for (size_t i = 1; i < 3; ++i)
	for (size_t j = 1; j < 3; ++j)
	    for (size_t k = 1; k < 3; ++k){
		Poly p;
		p.push_back(i);
		p.push_back(1);
		Poly q;
		q.push_back(j);
		q.push_back(1);
		Poly r;
		r.push_back(k);
		r.push_back(1);
		cout << "(" << p << ")"
		     << "(" << q << ")"
		     << "(" << r << ")"
		     << " &=& " << p * q * r << "\\\\" << endl;
	    }
    vector<Poly> q;
    q.push_back(Poly("101"));
    q.push_back(Poly("122"));
    q.push_back(Poly("112"));
    for (size_t i = 1; i < 3; ++i){
	Poly p;
	p.push_back(i);
	p.push_back(1);
	for (size_t j = 0; j < q.size(); ++j)
	    cout << "(" << q[j] << ")"
		 << "(" << p << ")"
		 << " &=& " << p * q[j] << "\\\\" << endl;
    }
}

void
factor()
{
    vector<Poly> p;
    p.push_back(Poly("1021"));
    p.push_back(Poly("1022"));
    p.push_back(Poly("1102"));
    p.push_back(Poly("1112"));
    p.push_back(Poly("1121"));
    p.push_back(Poly("1201"));
    p.push_back(Poly("1211"));
    p.push_back(Poly("1222"));
    p.push_back(Poly("10"));
    p.push_back(Poly("11"));
    p.push_back(Poly("12"));
    Poly q;
    q.push_back(1);
    for (size_t i = 0; i < p.size(); ++i){
	q = q * p[i];
	cout << p[i] << "\\\\" << endl;
    }
    cout << " = " << q << endl;
}

class Permutation : public vector<size_t>{
    
};

Permutation
operator * (const Permutation& g1, const Permutation& g2)
{
    Permutation ret;
    ret.resize(g1.size());
    for (size_t i = 0; i < g1.size(); ++i)
	ret[i] = g2[g1[i]];
    return ret;
}

void
buildCycles(const Permutation& g, vector<vector<size_t> >& rCycles)
{
    set<size_t> inCycle;
    for (size_t i = 0; i < g.size(); ++i){
	vector<size_t> cycle;
	for(size_t j = i; inCycle.count(j) == 0; j = g[j]){
	    cycle.push_back(j);
	    inCycle.insert(j);
	}
	
	if (cycle.size() > 1)
	    rCycles.push_back(cycle);
    }
}

ostream&
operator << (ostream& os, const Permutation& g)
{
    vector<vector<size_t> > cycles;
    buildCycles(g, cycles);
    for (size_t i = 0; i < cycles.size(); ++i){
	os << "(";
	for (size_t j = 0; j < cycles[i].size(); ++j)
	    os << cycles[i][j];
	os << ")";
    }
//    for (size_t j = 0; j < g.size(); ++j)
//	os << g[j];
    return os;
}

bool
even(const Permutation& g)
{
    vector<vector<size_t> > cycles;
    buildCycles(g, cycles);
    size_t numTrans = 0;
    for (size_t i = 0; i < cycles.size(); ++i)
	numTrans += cycles[i].size() + 1;
    return (numTrans % 2) == 0;
}

typedef set<Permutation> PermGroup;

PermGroup
generateA5()
{
    PermGroup A5;
    Permutation g;
    g.resize(5);
    for (size_t i = 0; i < 5; ++i)
	g[i] = i;
    for(; A5.count(g) == 0; next_permutation(g.begin(), g.end()))
	if (even(g))
	    A5.insert(g);
    return A5;
}

PermGroup
generate(const vector<Permutation>& v)
{
    PermGroup G;
    queue<Permutation> q;
    for (size_t i = 0; i < v.size(); ++i){
	G.insert(v[i]);
	q.push(v[i]);
    }
    while(q.size()){
	Permutation g = q.front();
	q.pop();
	for (PermGroup::const_iterator it = G.begin(); it != G.end(); ++it){
	    Permutation a = g * (*it); 
	    Permutation b = (*it) * g; 
	    if (G.insert(a).second)
		q.push(a);
	    if (G.insert(b).second)
		q.push(b);
	}
    }
    return G;
}

set<PermGroup>
addGenerator(const set<PermGroup>& groupSet, const PermGroup& G)
{
    set<PermGroup> ret;
    for (set<PermGroup>::const_iterator itS = groupSet.begin();
	 itS != groupSet.end(); ++itS)
    {
	for (PermGroup::const_iterator itG = G.begin(); itG != G.end(); ++itG){
	    vector<Permutation> generators;
	    generators.push_back(*itG);
	    generators.insert(generators.end(), itS->begin(), itS->end());
	    PermGroup H = generate(generators);
	    if (ret.insert(H).second)
		cout << H.size() << " ";
	}
    }
    return ret;
    
}

void
classifyA5()
{
    const PermGroup A5 = generateA5();
    set<PermGroup> singles;
    for (PermGroup::const_iterator it = A5.begin(); it != A5.end(); ++it){
	vector<Permutation> generators;
	generators.push_back(*it);
	singles.insert(generate(generators));
    }
    cout << "singles:" << endl << "    ";
    for (set<PermGroup>::const_iterator it = singles.begin(); 
	 it != singles.end(); ++it)
    {
	cout << it->size() << " ";
    }
    cout << endl;
    cout << "doubles:" << endl << "    ";
    set<PermGroup> doubles = addGenerator(singles, A5);
    cout << endl;
    cout << "triples:" << endl << "    ";
    set<PermGroup> triples = addGenerator(doubles, A5);
    cout << endl;
}



typedef Matrix<int> GL33El;

GL33El::Base
det(const GL33El& G)
{
    return 
	G(0, 0) * (G(1, 1) * G(2, 2) - G(1, 2) * G(2, 1)) -
	G(0, 1) * (G(1, 0) * G(2, 2) - G(2, 0) * G(1, 2)) +
	G(0, 2) * (G(1, 0) * G(2, 1) - G(1, 1) * G(2, 0));
}

// 012 
//0xxx
//1xxx
//2xxx

ostream& 
operator << (ostream& os, const GL33El& G)
{
    for (size_t i = 0; i < 3; ++i){
	for (size_t j = 0; j < 3; ++j)
	    os << G(i, j) << " ";
	os << endl;
    }
    return os;
}

void
normalize(GL33El& rG)
{
    for (size_t i = 0; i < 3; ++i)
	for (size_t j = 0; j < 3; ++j)
	    rG(i, j) = ((rG(i, j) % 3) + 3) % 3;
}

set<GL33El>
buildGL33()
{
    GL33El G(3, 3);
    set<GL33El> GL33;
    for (size_t n = 0; n < 19683; ++n){
	size_t m = n;
	for (size_t i = 0; i < 3; ++i)
	    for (size_t j = 0; j < 3; ++j){
		G(i, j) = m;
		m = m / 3;
	    }
	normalize(G);
	if (det(G) % 3 != 0)
	    GL33.insert(G);
     }
    cout << GL33.size() << endl; 
    return GL33;
}

GL33El
inv(const GL33El& m)
{
    GL33El md(m.rows(), 2 * m.cols());
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = 0; j < md.cols() / 2; ++j)
            md(i, j) = m(i, j);
    for (size_t i = 0; i < md.rows(); ++i)
        for (size_t j = md.cols() / 2; j < md.cols(); ++j)
            md(i, j) = 
                (i == j - md.cols() / 2) ? 1 : 0;
    for (size_t iter = 0; iter < min(m.rows(), m.cols()); ++iter){
        size_t ipiv = iter;
        for (; ipiv < md.rows() && (md(ipiv, iter) % 3) == 0; ++ipiv)
            ;
        assert(ipiv != md.rows());
        if (ipiv != iter)
            swapRows(md, ipiv, iter);
	int factor = md(iter, iter);
	for (size_t i = 0; i < md.cols(); ++i)
	    md(iter, i) = md(iter, i) * factor;
        for (size_t i = 0; i < md.rows(); ++i)
            if (i != iter)
                subRows(md, i, iter, md(i, iter));
    }
    
    GL33El minv(m.rows(), m.cols());
    for (size_t i = 0; i < minv.rows(); ++i)
        for (size_t j = 0; j < minv.cols(); ++j)
            minv(i, j) = md(i, j + md.cols() / 2);
    
//    assert(minv * m == Matrix<F>::eye(m.rows()));
    return minv;
}

void
findOrders()
{
    GL33El G(3, 3);
    G(0, 0) = 2; G(0, 1) = 1; G(0, 2) = 0;
    G(1, 0) = 0; G(1, 1) = 2; G(1, 2) = 1;
    G(2, 0) = 0; G(2, 1) = 0; G(2, 2) = 2;
    vector<GL33El> Gs;
    Gs.push_back(G);
    cout << "G=" << endl << G;
    for (size_t i = 1; i < 6; ++i){
	Gs.push_back(Gs[i - 1] * G);
	normalize(Gs.back());
	cout << "G^" << i + 1 << "=" << endl << Gs.back();
    }
}

void
disectGL33()
{
    
    const set<GL33El> GL33 = buildGL33();
    set<GL33El> todo = GL33;
    vector<set<GL33El> > classes;
    while (todo.size()){
	GL33El G = *todo.begin();
	classes.push_back(set<GL33El>());
	for (set<GL33El>::const_iterator it = GL33.begin(); it != GL33.end(); 
	     ++it)
	{
	    GL33El X = *it;
	    GL33El Y = inv(X);
	    GL33El Z = Y * X;
	    normalize(Y);
	    normalize(Z);
	    if (!(Z == GL33El::eye(3)))
		cout << Y << endl << X << endl << Z << endl << "xxx" << endl;
	    GL33El H = Y * G * X;
	    normalize(H);
	    if (classes.back().insert(H).second 
		&& H(1, 0)  == 0 && H(2, 0) == 0 
		&& H(2, 1)  == 0 && H(0, 2) == 0
		&& H(0, 1) < 2 && H(1, 2) < 2
		&& H(0, 0) <= H(1, 1)
		&& H(1, 1) <= H(2, 2))
	    {
/*		cout << H << "-----" << endl
		     << G << "-----" << endl
		     << X << "-----" << endl
		     << Y << "-----" << endl
		     << endl;
*/
	    }

	    todo.erase(H);
	}
	cout << classes.back().size() <<  " " << endl;
	for (set<GL33El>::const_iterator it = classes.back().begin();
	     it != classes.back().end(); ++it)
	{
	    GL33El X = *it;
	    if (it == classes.back().begin() ||
		(X(1, 0)  == 0 && X(2, 0) == 0 
		 && X(2, 1)  == 0 && X(0, 2) == 0
//		 && X(0, 0) <= X(1, 1)
//		 && X(1, 1) <= X(2, 2)
		 && ((X(0, 1) == 0 && X(1, 2) == 0) ||
		     (X(0, 1) == 1 && X(1, 2) == 0 && X(0, 0) == X(1, 1)) ||
		     (X(0, 1) == 0 && X(1, 2) == 1 && X(1, 1) == X(2, 2)) ||
		     (X(0, 1) == 1 && X(1, 2) == 1 && X(1, 1) == X(2, 2) && X(0, 0) == X(1, 1)))))
	    {
		cout << *it << endl;
	    }
	}

	cout << "=========================" << endl;
    }
    for (size_t i = 0; i < classes.size(); ++i)
	cout << classes[i].size() << " ";
    cout << endl;
}
int
main(int argc, char* argv[])
{
//    rootTable();
//    polyMath(argc, argv);
//    cubicPolys();
//    factor();
//    minpolyTable();
//    classifyA5();
    findOrders();
    disectGL33();
    return 0;
}
