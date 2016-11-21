#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "definitions.h"
#include "graph.h"
#include <tbb/tbb.h>

using namespace std;
#define EdgeSeparator SEPARATOR /**< Token used to separate strings and ints */


bool matchCompareUnique(Match lhs, Match rhs) {
    if (lhs.value == rhs.value) {
        if (lhs.lnode == rhs.lnode) {
            return (lhs.rnode == rhs.rnode);
        }
        return (lhs.lnode == rhs.lnode);
    } else {
        return lhs.value == rhs.value;
    }
}

bool matchCompareSort(Match lhs, Match rhs) {
    if (lhs.value == rhs.value) {
        if (lhs.lnode == rhs.lnode) {
            return (lhs.rnode > rhs.rnode);
        }
        return (lhs.lnode > rhs.lnode);
    } else {
        return lhs.value > rhs.value;
    }
}

Graph::Graph(bool directed) {
    _nodes = 0;
    _nEdges = NULL;
    _edges = NULL;
    __edges__ = NULL;

    _directed = directed;
}

Graph::Graph(const Graph &g) {

    _nodes = g._nodes;
    _nEdges = new int[_nodes];
    _edges = new int *[_nodes];

    int total_edges = 0; //Total edges is used to create a contiguous block of memory

    for (int i = 0; i < _nodes; i++) {
        _nEdges[i] = g._nEdges[i];
        total_edges += _nEdges[i];
    }

    __edges__ = new int[total_edges];

    total_edges = 0;
    for (int i = 0; i < _nodes; i++) {
        this->_edges[i] = __edges__ + total_edges; //Setting the pointers to the list of edges to the right place
        for (int j = 0; j < _nEdges[i]; j++) //Copying the edges
        {
            this->_edges[i][j] = g._edges[i][j];
        }
        total_edges += _nEdges[i];
    }

    _directed = g._directed;

}

Graph::Graph(const Graph &g, double s) //Create a graph by sampling from graph g by sampling probability s.
{

    for (int i = 0; i < g._nodes; i++) {
        _nodeIDtoint[i] = i;
        _inttonodeID[i] = i;
    }

    std::list<Edge> le;
    for (int i = 0; i < g._nodes; i++) {
        for (int j = 0; j < g._nEdges[i]; j++) //Copying the edges
        {
            if (g._edges[i][j] > i && ((double(rand()) / RAND_MAX) < s)) {
                Edge e;
                e.u = i;
                e.v = g._edges[i][j];
                le.push_back(e);
            }
        }
    }

    createGraph(g._nodes, le);

    le.clear();
}

Graph::Graph(const Graph &g, double s, double t) //Create a graph by sampling from graph g by sampling probability s.
{

    std::list<Edge> le;
    int nodesize = 0;
    for (int i = 0; i < g._nodes; i++) {
        if ((double(rand()) / RAND_MAX) < t) {
            _nodeIDtoint[i] = nodesize;
            _inttonodeID[nodesize] = i;
            nodesize += 1;
        }
    }
    int iID = 0;
    int jID = 0;

    cout << "\tgenerate subgraph" << endl;
    for (int i = 0; i < g._nodes; i++) {

        if (_nodeIDtoint.find(i) != _nodeIDtoint.end()) {
            for (int j = 0; j < g._nEdges[i]; j++) //Copying the edges
            {
                if (g._edges[i][j] > i && ((double(rand()) / RAND_MAX) < s) &&
                    _nodeIDtoint.find(g._edges[i][j]) != _nodeIDtoint.end()) {
                    Edge e;
                    e.u = _nodeIDtoint[i];
                    e.v = _nodeIDtoint[g._edges[i][j]];
                    le.push_back(e);
                }
            }
        }
    }
    cout << "\tV_count = " << _nodeIDtoint.size() << endl;
    createGraph(_nodeIDtoint.size(), le);

    le.clear();
}

inline std::string Graph::createMatchID(int i, int j) {
    std::stringstream sstm;
    sstm << i << "$" << j;
    return sstm.str();
}

void Graph::untokenizeMatchID(std::string pair, int &first, int &second) {
    std::string ES = EdgeSeparator;
    int index = pair.find(ES);
    first = atoi(pair.substr(0, index).c_str());
    if (index != -1)
        second = atoi(pair.substr(pair.find(ES) + ES.size()).c_str());
}

void Graph::createGraph(int nodes, list<Edge> &le) {
    initialize(nodes);
    //Calculating the number of neighbors of each node
    for (list<Edge>::iterator it = le.begin(); it != le.end(); ++it) {
        _nEdges[it->u]++;
        if (!_directed) {
            _nEdges[it->v]++;
        }
    }

    //Create a contiguous memory buffer for storing the edges
    long long size = 0;
    for (int i = 0; i < nodes; i++) {
        size += _nEdges[i];
    }
    __edges__ = new int[size];
    size = 0;
    for (int i = 0; i < nodes; i++) {
        _edges[i] = __edges__ + size;
        size += _nEdges[i];
    }

    //Read for the second time: writing edges: No need for try-catch block. File exists unless user deleted it...
    int *pos = new int[nodes]; //How many edges have been added for each node
    for (int i = 0; i < nodes; i++) {
        pos[i] = 0;
    }

    for (list<Edge>::iterator it = le.begin(); it != le.end(); ++it) //Writing the edges
    {
        int u = it->u, v = it->v;
        _edges[u][pos[u]++] = v;
        if (!_directed) {
            _edges[v][pos[v]++] = u;
        }
    }

    for (int i = 0; i < _nodes; i++) //Sort the neighbors to simplify search in the future
    {
        if (_nEdges[i] > 0) {
            random_shuffle(_edges[i], _edges[i] + _nEdges[i]);
            //sort( _edges[ i ], _edges[ i ] + _nEdges[ i ] );
        }
    }

    delete[] pos;
}

Graph::~Graph() {
    if (_nodes > 0) {
        delete[] __edges__;
        delete[] _edges;
        delete[] _nEdges;
    }
}

void Graph::initialize(int n) {
    _nodes = n;
    _nEdges = new int[n];
    _edges = new int *[n];

    for (int i = 0; i < n; i++) {
        _nEdges[i] = 0;
        _edges[i] = NULL;
    }
    _directed = false;
}

std::map<int, int> Graph::getnodeIDtoint() {
    return _nodeIDtoint;
}

std::map<int, int> Graph::getinttonodeID() {
    return _inttonodeID;
}

int Graph::getNEdges() {
    int ret = 0;
    for (int i = 0; i < _nodes; i++) {
        ret += _nEdges[i];
    }
    if (_directed) //If it is not a directed graph, then it means that each edge was counted twice
    {
        return ret;
    }
    return ret / 2;
}

void Graph::readGraph(string filename, bool is_shuffled) {
    int u, v, n;

    if (__edges__ != 0) //If the graph has been already constructed, avoid a memory leak
    {
        cout << __FILE__ << " " << __LINE__ << endl;
        throw ALLOC_EXCEPTION;
    }

    ifstream file;

    try {
        //Read for the first time: discover number of edges
        file.open(filename.c_str());
        file >> n; //Read the number of nodes
        initialize(n); //Initialize the graph with its correct size
        cout << "n = " << n << endl;
        if (is_shuffled) {
            while (file >> u >> v) {
                _nEdges[u - 1]++;
                if (!_directed) {
                    _nEdges[v - 1]++;
                }
            }
        } else {
            while (file >> u >> v) {

                _nEdges[u]++;
                if (!_directed) {
                    _nEdges[v]++;
                }
            }
        }


        file.close();
    }
    catch (...) {
        cout << "Problem reading the graph file: " << __FILE__ << " " << __LINE__ << endl;
        return;
    }

    //Create a contiguous memory buffer for storing the edges
    int size = 0;
    for (int i = 0; i < n; i++) {
        size += _nEdges[i];
    }
    __edges__ = new int[size];
    size = 0;
    for (int i = 0; i < n; i++) {
        _edges[i] = __edges__ + size;
        size += _nEdges[i];
    }

    //Read for the second time: writing edges: No need for try-catch block. File exists unless user deleted it...
    int *it = new int[n]; //How many edges have been added for each node
    for (int i = 0; i < n; i++) {
        it[i] = 0;
    }

    file.open(filename.c_str());
    file >> n;
    if (is_shuffled) {
        while (file >> u >> v) {
            u--;
            v--;
            _edges[u][it[u]++] = v;
            if (!_directed) //If it is not a directed graph, add the opposite edge too
            {
                _edges[v][it[v]++] = u;
            }
        }
    } else {
        while (file >> u >> v) {
            _edges[u][it[u]++] = v;
            if (!_directed) //If it is not a directed graph, add the opposite edge too
            {
                _edges[v][it[v]++] = u;
            }
        }
    }

    file.close();

    for (int i = 0; i < _nodes; i++) //Sort the neighbors to simplify search in the future
    {
        if (_nEdges[i] > 0)
            sort(_edges[i], _edges[i] + _nEdges[i]);
    }

    delete[] it;
}

int Graph::getNNodes() {
    return _nodes;
}

int Graph::getNNeighbors(int index) {
    if (index >= _nodes) {
        return 0;
    }
    return _nEdges[index];
}

int *Graph::getNeighbors(int index) {
    if (index >= _nodes) {
        return NULL;
    }
    return _edges[index];
}

bool Graph::isDirected() {
    return _directed;
}

bool Graph::edgeExists(int u, int v) //O( logN )
{
    if (_nEdges[u] == 0 or _nEdges[v] == 0) return false;
    int pos = binary_search(_nEdges[u], _edges[u], v);
    if (pos < 0 or pos >= _nodes) {
        cout << __FILE__ << " " << __LINE__ << endl;
        throw SEGFAULT;
    }
    return _edges[u][pos] == v;
}

void Graph::writeGraph(string output_file) {
    ofstream output(output_file.c_str());
    output << _nodes << endl; //Output number of nodes
    for (int i = 0; i < _nodes; i++) {
        for (int j = 0; j < _nEdges[i]; j++) {
            if (i > _edges[i][j]) //Prints each edge only once
                continue;
            output << i << " " << _edges[i][j] << endl;
        }
    }
    output.close();
}

double Graph::compareGraph(const Graph &g) //O( M )
{
    //Compare two graphs according to the existence (or not) of an edge in both graphs

    int total_size = 0, total_match = 0;
    int *it_this, *it_g;

    for (int u = 0; u < _nodes; u++) //O( M )
    {
        it_this = _edges[u];
        it_g = g._edges[u];

        total_size += (_nEdges[u] + g._nEdges[u]);

        //The neighbors of both graphs are sorted, so we can perform a matching in linear time
        while (it_this != _edges[u] + _nEdges[u] and it_g != g._edges[u] + g._nEdges[u]) {
            if (*it_this < *it_g) {
                it_this++;
            } else if (*it_this > *it_g) {
                it_g++;
            } else {
                it_this++;
                it_g++;
                total_match += 2;
            }
        }
    }

    if (total_size == 0)
        return 1.; //If there is no edge in both graphs, then it means that they match exactly, as they have the same number of neighbors

    return double(total_match) / double(total_size);
}

void Graph::display_time(const char *str) {
    time_t rawtime;
    time(&rawtime);
    cout << str << ": " << ctime(&rawtime);
}

std::string Graph::retrieveString(char *buf) {

    size_t len = 0;
    while (buf[len] != '\0') {
        len++;
    }
    return string(buf, len);
}

void Graph::showInterResults(vector<Match> &matches, Graph &g) {
    int correct = 0, wrong = 0;

    for (match_iterator it = matches.begin(); it != matches.end(); ++it) {
        if (_inttonodeID[it->lnode] == g._inttonodeID[it->rnode]) {
            correct++;
        } else {
            wrong++;
        }
    }
    cout << "Correct: " << correct << " Wrong: " << wrong << endl;
}

inline void
Graph::splitMarkVectors(RankedInactivePairs &ranked, MatchMultiset &inactive_pairs,
                        tbb::concurrent_unordered_set<int> &leftNodesM,
                        tbb::concurrent_unordered_set<int> &rightNodesM) {

    pair_iterator it_inctv;

    // A / M (just from begining) for detecting max MARKS
    int cur_markcount = 0;
    int r_indx = -1;
    for (it_inctv = inactive_pairs.begin(); it_inctv != inactive_pairs.end() && it_inctv->value >= 3; it_inctv++) {
        if (not_in_match(it_inctv->lnode, it_inctv->rnode)) {
            if (cur_markcount != it_inctv->value) {
                cur_markcount = it_inctv->value;
                r_indx++;
            }
            Match m;
            m.lnode = it_inctv->lnode;
            m.rnode = it_inctv->rnode;
            m.value = it_inctv->value;
            ranked[r_indx].push_back(m);
        }
    }

    // A / A(r >= 4)
    inactive_pairs.erase(inactive_pairs.begin(), it_inctv); //ToDo repair error

}

vector<Match> Graph::expandWhenStuckParallel(Graph &g, vector<Match> &seed_init, int r, int push_level, int push_size) {
    vector<Match> matches; // M

    // Structure that defines hashing and comparison operations for user's type.
    // PairTable global_score_map;
    // Localscore  local_score_map;
    PairTable score_map;

    MatchMultiset inactive_pairs; // marks for every pair mark count > r
    tbb::concurrent_unordered_set<int> leftNodesM; //G1(M)
    tbb::concurrent_unordered_set<int> rightNodesM; //G2(M)
    tbb::concurrent_unordered_set<string> seed_used; //Z
    tbb::concurrent_vector<Match> seed;

    unsigned int matchSize = 0;
    int push = 0; // for check push level ToDo what for is this?

    // M <- A_0 ps: A = A_0
    for (auto it_seed = seed_init.begin(); it_seed != seed_init.end(); ++it_seed) {
        Match m;
        m.lnode = it_seed->lnode;
        m.rnode = it_seed->rnode;
        m.value = 1;
        matches.push_back(m);
        leftNodesM.insert(it_seed->lnode);
        rightNodesM.insert(it_seed->rnode);
        seed.push_back(m);
    }

    // while |A| > 0 do
    int iter_num = 1;

    while (seed.size() > 0 && push < push_level) {
        cout << "############## iteration = " << iter_num++ << endl;
        cout << "seed size = " << seed.size() << endl;
        //for all pairs [i,j] of A do
        //ToDo begin 1 mapreduce
        typedef tbb::concurrent_vector<Match>::iterator concvec_it;
        tbb::parallel_for(tbb::blocked_range<concvec_it>(seed.begin(), seed.end()),
                          [&](tbb::blocked_range<concvec_it> &bl) {
                              auto mv = bl.begin();

                              int left_node = mv->lnode;
                              int right_node = mv->rnode;
                              string ID_i_j = createMatchID(left_node, right_node);
                              // add [i,j] to Z
                              seed_used.insert(ID_i_j);
                              // add one mark to all neighboring pairs of [i,j]
                              for (int i = 0; i < _nEdges[left_node]; i++) {
                                  for (int j = 0; j < g._nEdges[right_node]; j++) {
                                      if (not_in_matched(left_node, i, right_node, j)) {
                                          string ID_neighbor = createMatchID(_edges[left_node][i],
                                                                             g._edges[right_node][j]);
                                          PairTable::accessor a;
                                          if (score_map.insert(a, ID_neighbor)) {
                                              a->second = 1;
                                          } else { // update score
                                              a->second++;
                                          }
                                          a.release();
                                      }
                                  }
                              }
                          });
        //ToDo Reduce
        seed.clear();
        for (auto p_it = score_map.begin(); p_it != score_map.end(); p_it++) {
            Match m;
            untokenizeMatchID(p_it->first, m.lnode, m.rnode);
            m.value = p_it->second;
            inactive_pairs.insert(m);
        }
        score_map.clear();

        //ToDo end 1 MapReduce
        // A <- None

        cout << "MapReduce_1 end" << endl;

        //ToDo deadlock is after somewhere there

        bool stopFlag = true; //ToDo what is this about
        // while there exists an unmatched pair with score at least r+1
        while (inactive_pairs.size() > 0 && (inactive_pairs.begin()->value > r) && matches.size() < 200000000) {

            //get ranked r>=4
            RankedInactivePairs ranked;
            splitMarkVectors(ranked, inactive_pairs, leftNodesM, rightNodesM);
            vector<size_t> vec;
            for (size_t i = 0; i < ranked.size(); ++i) {
                vec.push_back(i);
            }

            tbb::concurrent_vector<Match> local_matches;
            //ToDo begin 2 MapReduce
            tbb::parallel_for_each(vec.begin(), vec.end(),
                                   [&local_matches, &leftNodesM, &rightNodesM, &ranked, &g, this](size_t r_indx) {
                                       auto vecA = ranked[r_indx];
                                       for (auto it_inctv = vecA.begin(); it_inctv != vecA.end(); it_inctv++) {
                                           int left_node = it_inctv->lnode;
                                           int right_node = it_inctv->rnode;
                                           if (not_in_match(left_node, right_node)) {  //EXC_BAD_ACCESS
                                               if (abs(_nEdges[left_node] - g._nEdges[right_node]) == 0) {
                                                   Match mtchd;
                                                   mtchd.lnode = left_node;
                                                   mtchd.rnode = right_node;
                                                   mtchd.value = 0;
                                                   local_matches.push_back(mtchd);
                                                   leftNodesM.insert(left_node);
                                                   rightNodesM.insert(right_node);
                                               }
                                           }
                                       }
                                   });

            // among the pairs with the highest score
            auto it_inctv = inactive_pairs.begin();
            // A / M for detecting max MARKS
            while (inactive_pairs.size() > 0) {
                if (not_in_match(it_inctv->lnode, it_inctv->rnode)) {
                    break; // TRANSFER TO
                } else if (it_inctv != inactive_pairs.end()) {
                    inactive_pairs.erase(it_inctv++);
                }
            }
            // check if |A| > 0
            if (inactive_pairs.size() == 0) {
                break;
            }

            for (int marks_count = it_inctv->value; marks_count >= 2; marks_count--) {
                // select the unmatched pair [i,j] with the minimum |d_1,i - d_2,j| (degree diff)
                cout << inactive_pairs.begin()->lnode << " " << inactive_pairs.begin()->rnode << " "
                     << inactive_pairs.begin()->value << endl;
                for (it_inctv = inactive_pairs.begin();
                     it_inctv != inactive_pairs.end() && it_inctv->value == marks_count; it_inctv++) {
                    //ToDo EXC_BAD_ACCESS
                    // select the unmatched pair [i,j]
                    if (not_in_match(it_inctv->lnode, it_inctv->rnode)) {
                        // find min |d_1,i - d_2,j|
                        if (abs(_nEdges[it_inctv->lnode] - g._nEdges[it_inctv->rnode]) == 0) {

                            int left_node = it_inctv->lnode;
                            int right_node = it_inctv->rnode;
                            Match m;
                            m.lnode = left_node;
                            m.rnode = right_node;
                            m.value = 1;
                            local_matches.push_back(m);
                            leftNodesM.insert(left_node);
                            rightNodesM.insert(right_node);
                        }
                    }
                    inactive_pairs.erase(it_inctv);
                }
            }
            //ToDo end 2 MapReduce
            cout << "MapReduce_2 end" << endl;

            //ToDo start 3 MapReduce
//            tbb::parallel_for_each(local_matches.begin(), local_matches.end(),
//                                   [&](Match &m) {
            tbb::parallel_for(tbb::blocked_range<concvec_it>(local_matches.begin(), local_matches.end()),
                              [&](tbb::blocked_range<concvec_it> &bl) {
                                        auto m = bl.begin();
                                       int left_node = m->lnode;
                                       int right_node = m->rnode;
                                       string ID_not_active = createMatchID(left_node, right_node);
                                       // if [i,j] not in Z
                                       if (seed_used.find(ID_not_active) == seed_used.end()) { //ToDo change to contains
                                           // add [i,j] to Z
                                           seed_used.insert(ID_not_active);
                                           // add one marks to all of its neighbouring pairs
                                           for (int i = 0; i < _nEdges[left_node]; i++) {
                                               for (int j = 0; j < g._nEdges[right_node]; j++) {
                                                   if (not_in_matched(left_node, i, right_node, j)) {
                                                       string ID_neighbor = createMatchID(_edges[left_node][i],
                                                                                          g._edges[right_node][j]);
                                                       PairTable::accessor a;
                                                       if (score_map.insert(a, ID_neighbor)) {
                                                           a->second = 1;
                                                       } else { // update score
                                                           a->second++;
                                                       }
                                                       a.release();
                                                   }
                                               }
                                           }
                                       }
                                   });
            //ToDo Reduce

            for (auto p_it = score_map.begin(); p_it != score_map.end(); p_it++) {
                Match m;
                untokenizeMatchID(p_it->first, m.lnode, m.rnode);
                m.value = p_it->second;
                inactive_pairs.insert(m);
            }
            score_map.clear();
            //ToDo end 3 MapReduce
            cout << "MapReduce_3 end" << endl;

            //ToDo cp local_matches to matches
            for (auto m_it = local_matches.begin(); m_it != local_matches.end(); m_it++) {
                Match m;
                m.lnode = m_it->lnode;
                m.rnode = m_it->rnode;
                matches.push_back(m);
            }

        }
        cout << "Finish with inactive_pairs" << endl;

        int seedCnt = 0;
        cout << "++++++" << endl;

        push++;
        //ToDo begin 4 MapReduce
        // A <- all neighbors of M [i,j] not in Z, i,j not in V_1,V_2(M)
        cout << "MapReduce 4 start" << endl;
        tbb::parallel_for(tbb::blocked_range<match_iterator>(matches.begin(), matches.end()),
                          [&, this](tbb::blocked_range<match_iterator> &bl) {
                              auto it = bl.begin();
                              int left_node = it->lnode;
                              int right_node = it->rnode;
                              // all neighbors of M
                              for (int i = 0; i < _nEdges[left_node]; i++) {
                                  for (int j = 0; j < g._nEdges[right_node]; j++) {
                                      // i,j not in V_1,V_2(M)
                                      if (not_in_matched(left_node, i, right_node, j)) {
                                          string ID = createMatchID(_edges[left_node][i], g._edges[right_node][j]);
                                          // [i,j] not in Z
                                          if (seed_used.find(ID) == seed_used.end()) { //ToDo && seedCnt < push_size
                                              Match m;
                                              m.lnode = _edges[left_node][i];
                                              m.rnode = g._edges[right_node][j];
                                              m.value = 1;
                                              seed.push_back(m); //ToDo SIGABRT
                                              seedCnt++;
                                          }
                                      }

                                  }
                              }
                          });
        //ToDo end 4 MapReduce

        cout << "Extended seed size (MapReduce4 end): " << seed.size() << endl;
//            if (push == push_level - 1 or seed.size() == 0) { //ToDo decide about this one
//
//                cout << "size inactive_pairs with value = 1: " << inactive_pairs.size() << endl;
//                for (auto it = inactive_pairs.begin(); it != inactive_pairs.end(); ++it) { //ToDo EXC_BAD_ACCESS
//                    if (it->value == 1) {
//                        inactive_pairs.erase(it);
//                    }
//                }
//                cout << "size: " << inactive_pairs.size() << endl;
//                inactive_pairs.clear();
//                cout << "Finish with inactive_pairs" << endl;
//            }
    }

    inactive_pairs.clear();
    seed.clear();
    leftNodesM.clear();
    rightNodesM.clear();
    return matches;
}