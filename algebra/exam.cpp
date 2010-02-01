#include <iostream>
#include <fstream>
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

int
main()
{
    rootTable();
    return 0;
}
