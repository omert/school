#include <stdlib.h>
#include <iostream>
#include <set>
#include <vector>
#include <map>

using namespace std;

typedef vector<size_t> EdgeSet;
typedef vector<EdgeSet> Graph;

void
populate(Graph& rG, size_t n)
{
    rG.resize(n);
    for (size_t i = 0; i < n; ++i)
	for (size_t j = 0; j < n; ++j)
	    if (i != j)
		rG[i].push_back(j);
}

//map<size_t, size_t>
size_t
walkSim(size_t n)
{
    Graph G;
    populate(G, n);

    size_t step = 0;
    size_t pos = 0;
    map<size_t, size_t> visits;
    size_t lastVisit = 0;
//    cout << pos;
    for(; G[pos].size(); ++step){
	size_t edgePos = rand() % G[pos].size();
	size_t nextPos = G[pos][edgePos];
	for (size_t k = 0; k < G[nextPos].size(); ++k)
	    if (G[nextPos][k] == pos){
		G[nextPos][k] = G[nextPos].back();
		G[nextPos].pop_back();
	    }
	G[pos][edgePos] = G[pos].back();
	G[pos].pop_back();
	pos = nextPos;
	if (nextPos == 0){
	    visits[visits.size()] = step - lastVisit;
	    lastVisit = step;
	}
	    
//	cout << " -> " << pos;
    }
    return step;
//    return visits;
//    cout << endl << step << " steps" << endl;
/*
    map<size_t, size_t> hist;
    for (size_t i = 0; i < n; ++i)
	hist[G[i].size()]++;
    return hist;
/*
    size_t numEmpty = 0;
    for (size_t i = 0; i < n; ++i)
	if (G[i].empty())
	    ++numEmpty;
    return numEmpty;
*/
}

int
main(int argc, char* argv[])
{
    if (argc < 3){
	cout << "2 parameter needed" << endl;
	exit(-1);
    }
    srand(atoi(argv[2]));
    size_t n = atoi(argv[1]);
    size_t sumSteps = 0;
    size_t iters = 10000;
    map<size_t, size_t> hist;
    for (size_t i = 0; i < iters; ++i){
	sumSteps += walkSim(n);
//	map<size_t, size_t> h = walkSim(n);
//	for(size_t j = 0; j < h.size(); ++j)
//	    hist[j] += h[j];

    }
//    for(size_t j = 0; j < hist.size(); ++j)
//	cout << j << ": " << 1.0 * hist[j] / iters  << endl;
    
    cout << (1.0 * n * (n - 1) / 2 - 1.0 * sumSteps / iters) / n << endl;
	
    return 0;
}
