#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <map>

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
    bool printed = false;
    for (size_t i = 0; i < p.size(); ++i)
	if (p[i] != 0){
	    if (printed)
		os << " + ";
	    printed = true;
	    if (i == 0)
		os << p[i];
	    else{
		if (p[i] != 1)
		    os << p[i];
		if (i == 1)
		    os << "x";
		else
		    os << "x^" << i;
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
    
}

int
main(int argc, char* argv[])
{
//    rootTable();
    polyMath(argc, argv);
    cubicPolys();
    return 0;
}
