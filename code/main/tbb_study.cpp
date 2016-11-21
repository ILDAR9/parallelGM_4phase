#include <tbb/tbb.h>
#include <iostream>
#include <cmath>
#include <set>
#include <map>

using namespace std;

void test1(void) {
    tbb::flow::graph g;
    tbb::flow::continue_node<tbb::flow::continue_msg>
            h(g, [](const tbb::flow::continue_msg &) { cout << "Hello "; });

    tbb::flow::continue_node<tbb::flow::continue_msg>
            w(g, [](const tbb::flow::continue_msg &) { cout << "Hello "; });

    tbb::flow::make_edge(h, w);

    h.try_put(tbb::flow::continue_msg());

    g.wait_for_all();
}

void test2(void) {
    int size = 20000000;
    double *input = new double[size];
    double *output = new double[size];


    for (int i = 0; i < size; ++i) {
        input[i] = i;
    }

    for (int i = 0; i < size; ++i) {
        output[i] = sqrt(sin(input[i]) * sin(input[i]) + cos(input[i]) * cos(input[i]));
    }
    cout << output[size - 1];

}

void test3(void) {
    int size = 20000000;
    double *input = new double[size];
    double *output = new double[size];


    for (int i = 0; i < size; ++i) {
        input[i] = i;
    }

    tbb::parallel_for(0, size, 1, [=](int i) {
        output[i] = sqrt(sin(input[i]) * sin(input[i]) + cos(input[i]) * cos(input[i]));
    });

}

void test4(void) {
    struct MyHashCompare {
        static size_t hash(const string &x) {
            size_t h = 0;
            for (const char *s = x.c_str(); *s; ++s)
                h = (h * 17) ^ *s;
            return h;
        }

        //! True if strings are equal
        static bool equal(const string &x, const string &y) {
            return x == y;
        }
    };
    typedef tbb::concurrent_hash_map<string, int, MyHashCompare> StringTable;
    StringTable table;

    vector<string> vec;
    int N = 1000;
    for (int i = 0; i < N; ++i) {
        vec.push_back(to_string(i));
    }

    tbb::parallel_do(vec.begin(), vec.end(), [&table](string &val) {
        StringTable::accessor a;
        table.insert(a, val);

        if (stoi(val) == 100) {
            a->second += 200;
        } else {
            a->second += 3;
        }
        a.release();
    });

    StringTable::accessor a;
    table.find(a, "100");


    cout << a->first << " " << a->second << endl;
    table.erase(a);
    a.release();
    cout << table.find(a, "200") << endl;
    cout << a->first << " " << a->second << endl;
    a.release();

}

void test5() {
    multiset<int> a;
    for (int i = 0; i < 10; ++i) {
        a.insert(i);
    }

    typedef multiset<int>::iterator mit;
    mit it;

    for (it = a.begin(); it != a.end() && *it < 0; it++);

    a.erase(a.begin(), it);

    for (it = a.begin(); it != a.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;

}

typedef map<int, vector<int>> m_vec;

void test6(m_vec &m) {
    m[2].push_back(23);
}

void test7() {
    typedef vector<int> v_t;
    v_t v;
    auto cp_match = [&v]() {
        v_t vc(v);

        return &vc;
    };

    v.push_back(23);

    v_t *t;
    t = cp_match();
    cout << *(*t).begin() << endl;

}

void test8() {
    struct Mo {
        int val;
    };

    struct classcomp {
        bool operator()(const Mo &lhs, const Mo &rhs) const { return lhs.val < rhs.val; }
    };

    typedef set<Mo, classcomp> stype;
    stype s;
    for (int i = 5; i >= 0; i--) {
        Mo m;
        m.val = i;
        s.insert(m);
    }
    for (stype::iterator it = s.begin(); it != s.end(); it++) {
        cout << it->val << " ";
        if (it->val == 2) {
            //ToDO
        }
    }
    cout << endl;

    for (auto it = s.begin(); it != s.end(); it++) {
        cout << it->val << " ";
    }

    cout << endl;
}

void test9_1() {
    std::map<int, int> m;
    auto size = 1000000;
    for (int i = 0; i < size; ++i) {
        m[i % 10000]  = i;
    }
    cout << m.size();

}

void test9_2() {
    typedef tbb::concurrent_hash_map<int, int> PairTable;
    PairTable score_map;

    auto size = 1000000;
    tbb::parallel_for(0, size, 1, [&](int i) {
        PairTable::accessor a;
        score_map.insert(a, i % 10000);
        a->second = 1;

        a.release();
    });
    cout << score_map.size();
}


int main(void) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//    m_vec m;
    test9_1();
//    cout << *m[3].begin() << endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << std::endl;

    return 0;
}
