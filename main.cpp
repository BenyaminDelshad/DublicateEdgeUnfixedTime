// Written By Benyamin Delshad Mamaghani - 2019
// We use sorted Temporal Edge Format (node1,node2, Time, landa = 1) in this code with time interval [0, oo]
// we have to dublicate one edge with our arbitrary time to improve diameter property.
// ignoring infinity diameter!
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

using namespace std;

const long long MAXN = 2000 + 10;
const long long MAXM = 100000 + 10;
const long long oo = 1000ll * 1000ll * 1000ll * 1000ll * 1000ll * 1000ll + 10;

class TemporalEdge{
public:
    long long node1, node2, t, landa;
    TemporalEdge() {
        node1 = node2 = t = landa = 0;
    }
    TemporalEdge(long long a, long long b, long long c, long long d) {
        node1 = a; node2 = b; t = c; landa = d;
    }
};


//string input_file_name = "../transportation/grenoble/grenoble-sorted.txt";
//string output_file_name = "../transportation/grenoble-output.txt";
string input_file_name = "../good_input_1.txt";
string output_file_name = "../good_output_1.txt";
//string input_file_name = "../good_input_2.txt";
//string output_file_name = "../good_output_2.txt";

vector<TemporalEdge>G; // Sorted Edge Stream
vector<TemporalEdge>Candidate_method2;
vector<TemporalEdge>Candidate_method3;
vector<TemporalEdge>Candidate_method4;

int par[MAXN][MAXN]; // ID OF EDGE!
int furthest_node[MAXN];
map<pair<int, int>, TemporalEdge>last_used_edge;
bool is_node[MAXN];
long long last_time[MAXN][MAXN]; // method 4 material
//vector<long long> input_e[MAXN]; // list of edges that go to this node
//vector<long long>output_e[MAXN]; // list of edges that go somewhere else from this node
long long last_start_time[MAXN]; // max(input of this node plus landa)
long long t_alpha, t_omega; // time long interval
long long m = 0; // number of edges
long long n = 1; // number of vertices
long long Original_diam;

ifstream fin(input_file_name);
ofstream fout(output_file_name);

inline bool operator< (const TemporalEdge& a, const TemporalEdge& b){ return a.t < b.t;};
inline bool operator== (const TemporalEdge& a, const TemporalEdge& b){ return (a.node1 == b.node1 && a.node2 == b.node2 && a.t == b.t && a.landa == b.landa);};
void TemporalEdge_eq (TemporalEdge& a, TemporalEdge& b){a.node1 = b.node1; a.node2 = b.node2; a.t = b.t; a.landa = b.landa;};

vector<long long> foremost_path(vector<TemporalEdge> &, long long, long long = n, long long = t_alpha, long long = t_omega);
long long calculate_diam(vector<TemporalEdge> &, long long = n, long long = t_alpha, long long = t_omega);
vector<long long> origin_foremost_path(vector<TemporalEdge> &, long long, long long = n, long long = t_alpha, long long = t_omega);
long long origin_calculate_diam(vector<TemporalEdge> &, long long = n, long long = t_alpha, long long = t_omega);

void input();
void initiate();
void method1(); // dublicate each edge with new time max of inputs plus landa.
void method2();// considering longest path from each node and dublicate last edge of each one with par plus landa
void method3(); // like method 2 but dublicate all edges of each longest path
void method4(); // like method 3 but use at most one edge for each pair of node(last one)

int main() {
    initiate();
    input();
    Original_diam = origin_calculate_diam(G, n, 0, oo);
    cout << Original_diam << endl;
    fout << "\nOriginal Diam: " << Original_diam << endl;
    //method1();
    method2();
    //method3();
    method4();
    return 0;
}
void method1 () { // dublicate each edge with new time max of inputs plus landa.
    for(int i = 0; i < n;++i)
        last_start_time[i] = 0;
    for(int i = 0; i < m; ++i) {
        last_start_time[G[i].node2] = max(last_start_time[G[i].node2], G[i].t + G[i].landa);
    }
    vector<TemporalEdge>G2;
    for(int i = 0; i < m; ++i) {
        G2.push_back(TemporalEdge(G[i].node1, G[i].node2, G[i].t, G[i].landa));
    }
    long long cur_diam = Original_diam;
    TemporalEdge choice(G2[0].node1, G2[0].node2, G2[0].t, G2[0].landa);
    set<pair<int, int> > has_been_checked;
    for(int i = 0; i < m; ++i) {
        if(has_been_checked.find(make_pair(G[i].node1, G[i].node2)) != has_been_checked.end())
            continue;
        has_been_checked.insert(make_pair(G[i].node1, G[i].node2));
        TemporalEdge candidate(G[i].node1, G[i].node2, last_start_time[G[i].node1], 1);
        G2.push_back(candidate);
        sort(G2.begin(), G2.end());
        long long tmp_diam = calculate_diam(G2);
        if(tmp_diam < cur_diam) {
            cur_diam = tmp_diam;
            choice = candidate;
        }
        for(int j = 0; j < m + 1; ++j) {
            if(candidate == G2[j]) {
                G2.erase(G2.begin() + j);
                break;
            }
        }
    }
    // now we have our choice and new diameter
    fout << "\n\nMethod1 (dublicate each edge with new time max of inputs plus landa.): \nNew edge: ";
    fout << choice.node1 << " " << choice.node2 << " " << choice.t << " " << choice.landa << endl;
    fout << "New Diameter: " << cur_diam << endl;
    return;
}
void method2() { // considering longest path from each node and dublicate last edge of each one with par plus landa
    for(int i = 0; i < n; ++i) {
        if(!is_node[i]) continue;
        int des = furthest_node[i];
        TemporalEdge tmp(0,0,0,0);
        TemporalEdge_eq(tmp, last_used_edge[make_pair(i, des)]);
        //cerr << tmp.node1 << ' ' << tmp.node2 << ' ' << tmp.t << ' ' << tmp.landa << endl;
        long long tmp_t = 0;
        if(par[i][tmp.node1] == -1) tmp_t = 0;
        else tmp_t = G[par[i][tmp.node1]].t + G[par[i][tmp.node1]].landa;
        Candidate_method2.push_back(TemporalEdge(tmp.node1, tmp.node2, tmp_t, 1));
    }

    vector<TemporalEdge>G2;
    for(int i = 0; i < m; ++i) {
        G2.push_back(TemporalEdge(G[i].node1, G[i].node2, G[i].t, G[i].landa));
    }
    long long cur_diam = Original_diam;
    TemporalEdge choice(Candidate_method2[0].node1, Candidate_method2[0].node2, Candidate_method2[0].t, Candidate_method2[0].landa);
    for(int i = 0; i < int(Candidate_method2.size()); ++i) {
        cerr << Candidate_method2[i].node1 << ' ' << Candidate_method2[i].node2 << ' ' << Candidate_method2[i].t << endl;
        TemporalEdge candidate(Candidate_method2[i].node1, Candidate_method2[i].node2, Candidate_method2[i].t, Candidate_method2[i].landa);
        G2.push_back(candidate);
        sort(G2.begin(), G2.end());
        long long tmp_diam = calculate_diam(G2);
        if(tmp_diam < cur_diam) {
            cur_diam = tmp_diam;
            choice = candidate;
        }
        for(int j = 0; j < m + 1; ++j) {
            if(candidate == G2[j]) {
                G2.erase(G2.begin() + j);
                break;
            }
        }
    }
    // now we have our choice and new diameter
    fout << "\n\nMethod2 (considering longest path from each node and dublicate last edge of each one with par plus landa.): \nNew edge: ";
    fout << choice.node1 << " " << choice.node2 << " " << choice.t << " " << choice.landa << endl;
    fout << "New Diameter: " << cur_diam << endl;
    return;
}
void method3() { // like method 2 but dublicate all edges of each longest path
    for(int i = 0; i < n; ++i) {
        if(!is_node[i]) continue;
        int des = furthest_node[i];
        int cur = des;
        //cerr << i << ' ' << cur << endl;
        while(i != cur) {
            TemporalEdge tmp(0,0,0,0);
            TemporalEdge_eq(tmp, last_used_edge[make_pair(i, cur)]);
            cur = G[par[i][cur]].node1;
            //cerr << tmp.node1 << ' ' << tmp.node2 << ' ' << tmp.t << ' ' << tmp.landa << endl;
            long long tmp_t = 0;
            if (par[i][tmp.node1] == -1) tmp_t = 0;
            else tmp_t = G[par[i][tmp.node1]].t + G[par[i][tmp.node1]].landa;
            Candidate_method3.push_back(TemporalEdge(tmp.node1, tmp.node2, tmp_t, 1));
        }
    }

    vector<TemporalEdge>G2;
    for(int i = 0; i < m; ++i) {
        G2.push_back(TemporalEdge(G[i].node1, G[i].node2, G[i].t, G[i].landa));
    }
    long long cur_diam = Original_diam;
    TemporalEdge choice(Candidate_method3[0].node1, Candidate_method3[0].node2, Candidate_method3[0].t, Candidate_method3[0].landa);
    cerr << "Method 3: \n\n";
    for(int i = 0; i < int(Candidate_method3.size()); ++i) {
        cerr << Candidate_method3[i].node1 << ' ' << Candidate_method3[i].node2 << ' ' << Candidate_method3[i].t << endl;
        TemporalEdge candidate(Candidate_method3[i].node1, Candidate_method3[i].node2, Candidate_method3[i].t, Candidate_method3[i].landa);
        G2.push_back(candidate);
        sort(G2.begin(), G2.end());
        long long tmp_diam = calculate_diam(G2);
        if(tmp_diam < cur_diam) {
            cur_diam = tmp_diam;
            choice = candidate;
        }
        for(int j = 0; j < m + 1; ++j) {
            if(candidate == G2[j]) {
                G2.erase(G2.begin() + j);
                break;
            }
        }
    }
    // now we have our choice and new diameter
    fout << "\n\nMethod3 (like method 2 but dublicate all edges of each longest path.): \nNew edge: ";
    fout << choice.node1 << " " << choice.node2 << " " << choice.t << " " << choice.landa << endl;
    fout << "New Diameter: " << cur_diam << endl;
    return;
}
void method4(){ // like method 3 but use at most one edge for each pair of node(last one)
    for(int i = 0; i < n; ++i) {
        if(!is_node[i]) continue;
        int des = furthest_node[i];
        int cur = des;
        //cerr << i << ' ' << cur << endl;
        while(i != cur) {
            TemporalEdge tmp(0,0,0,0);
            TemporalEdge_eq(tmp, last_used_edge[make_pair(i, cur)]);
            cur = G[par[i][cur]].node1;
            //cerr << tmp.node1 << ' ' << tmp.node2 << ' ' << tmp.t << ' ' << tmp.landa << endl;
            long long tmp_t = 0;
            if (par[i][tmp.node1] == -1) tmp_t = 0;
            else tmp_t = G[par[i][tmp.node1]].t + G[par[i][tmp.node1]].landa;
            Candidate_method4.push_back(TemporalEdge(tmp.node1, tmp.node2, tmp_t, 1));
        }
    }
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n ; ++j)
            last_time[i][j] = -1;
    for(int i = 0; i < int(Candidate_method4.size()); ++i) {
        last_time[Candidate_method4[i].node1][Candidate_method4[i].node2] = max(Candidate_method4[i].t, last_time[Candidate_method4[i].node1][Candidate_method4[i].node2]);
    }
    Candidate_method4.clear();
    for(int i = 0; i < n; ++i) {
        if(!is_node[i]) continue;
        for (int j = 0; j < n; ++j) {
            if(!is_node[j]) continue;
            if(last_time[i][j] == -1) continue;
            Candidate_method4.push_back(TemporalEdge(i,j,last_time[i][j], 1));
        }
    }
    vector<TemporalEdge>G2;
    for(int i = 0; i < m; ++i) {
        G2.push_back(TemporalEdge(G[i].node1, G[i].node2, G[i].t, G[i].landa));
    }
    long long cur_diam = Original_diam;
    TemporalEdge choice(Candidate_method4[0].node1, Candidate_method4[0].node2, Candidate_method4[0].t, Candidate_method4[0].landa);
    cerr << "Method 4: \n\n";
    for(int i = 0; i < int(Candidate_method4.size()); ++i) {
        cerr << Candidate_method4[i].node1 << ' ' << Candidate_method4[i].node2 << ' ' << Candidate_method4[i].t << endl;
        TemporalEdge candidate(Candidate_method4[i].node1, Candidate_method4[i].node2, Candidate_method4[i].t, Candidate_method4[i].landa);
        G2.push_back(candidate);
        sort(G2.begin(), G2.end());
        long long tmp_diam = calculate_diam(G2);
        if(tmp_diam < cur_diam) {
            cur_diam = tmp_diam;
            choice = candidate;
        }
        for(int j = 0; j < m + 1; ++j) {
            if(candidate == G2[j]) {
                G2.erase(G2.begin() + j);
                break;
            }
        }
    }
    // now we have our choice and new diameter
    fout << "\n\nMethod4 (like method 3 but use at most one edge for each pair of node(last one).): \nNew edge: ";
    fout << choice.node1 << " " << choice.node2 << " " << choice.t << " " << choice.landa << endl;
    fout << "New Diameter: " << cur_diam << endl;
    return;
}
void initiate() {
    t_alpha = 0;
    t_omega = oo;
    for(int i = 0; i < MAXN; ++i)
        is_node[i] = false;
    G.clear();
    return;
}
void input() {
    initiate();
    long long tmp1, tmp2, tmp_time, tmp_landa = 1;
    while(!fin.eof()) {
        fin >> tmp1 >> tmp2 >> tmp_time;
        is_node[tmp1] = is_node[tmp2] = true;
        G.push_back(TemporalEdge(tmp1, tmp2, tmp_time, tmp_landa));
        if(tmp1 + 1 > n) n = tmp1 + 1;
        if(tmp2 + 1 > n) n = tmp2 + 1;
    }
    m = G.size();
    return;
}
vector<long long> foremost_path(vector<TemporalEdge> &G2, long long source2, long long n2, long long t_alpha2, long long t_omega2){
    vector<long long>t;
    //for(long long i = 0; i < n2; ++i)
    // initialize t[v] for all v in V(G)
    for(long long i = 0; i < n2; ++i)
        t.push_back(oo);
    t[source2] = t_alpha2;
    // sorting edges
    // Edges sorted
    // body
    long long m2 = G2.size();
    for(long long i = 0; i < m2; ++i) {
        if(G2[i].t + G2[i].landa <= t_omega2 && G2[i].t >= t[G2[i].node1]) {
            if(G2[i].t + G2[i].landa < t[G2[i].node2]) {
                t[G2[i].node2] = G2[i].t + G2[i].landa;

            }
        }
        else if(G2[i].t >= t_omega2) {
            break;
        }
    }

    // t is OK!
    return t;
}

long long calculate_diam(vector<TemporalEdge> &G2, long long n2 , long long t_alpha2 , long long t_omega2) {
    //cout << n2 << endl;
    vector<long long>e;
    for(int i = 0; i < n2; ++i)
        e.push_back(0);
    long long diam = 0;
    int number_of_bad_pairs = 0;
    //cout << t_alpha2 << ' ' << t_omega2 << endl;
    for(int i = 0; i < n2; ++i) {
        if(!is_node[i]) continue;
        vector<long long>t = foremost_path(G2, i, n2, t_alpha2, t_omega2);
        //cout << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << ' ' << t[4] << ' ' << t[5] << endl;
        for(int j = 0; j < n2; ++j) {
            if (is_node[j] and t[j] < oo) {
                e[i] = max(e[i], t[j]);
                //if(i == 1056 and j == 771) cout << t[j] << endl; // OK!
                if (t[j] == oo) number_of_bad_pairs++;
            }
        }

        //cout << e[i] << endl;
    }

    //cout << number_of_bad_pairs << endl;
    //int number_nodes=0;
    //for(int j = 0; j < n; ++j){
    //    if(is_node[j]) number_nodes++;
    //}
    //cout << number_nodes << ' ' << number_of_bad_pairs << endl; // 1547, 1107093 out of 2393209 (46%)


    for(int i = 0; i < n2; ++i) {
        //if(e[i] >= oo) cout << i << endl;
        diam = max(diam, e[i]);
    }
    return diam;
}

vector<long long> origin_foremost_path(vector<TemporalEdge> &G2, long long source2, long long n2, long long t_alpha2, long long t_omega2){
    vector<long long>t;
    //for(long long i = 0; i < n2; ++i)
    // initialize t[v] for all v in V(G)
    for(long long i = 0; i < n2; ++i)
        t.push_back(oo);
    t[source2] = t_alpha2;
    // sorting edges
    // Edges sorted
    // body
    long long m2 = G2.size();
    for(long long i = 0; i < m2; ++i) {
        if(G2[i].t + G2[i].landa <= t_omega2 && G2[i].t >= t[G2[i].node1]) {
            if(G2[i].t + G2[i].landa < t[G2[i].node2]) {
                t[G2[i].node2] = G2[i].t + G2[i].landa;

                last_used_edge[make_pair(source2,G2[i].node2)] = TemporalEdge(G2[i].node1, G2[i].node2, G2[i].t, G2[i].landa);
                par[source2][G2[i].node2] = i;
            }
        }
        else if(G2[i].t >= t_omega2) {
            break;
        }
    }

    // t is OK!
    return t;
}

long long origin_calculate_diam(vector<TemporalEdge> &G2, long long n2 , long long t_alpha2 , long long t_omega2) {
    for(int i = 0 ; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            par[i][j] = -1;
        }
    }
    //cout << n2 << endl;
    vector<long long>e;
    for(int i = 0; i < n2; ++i)
        e.push_back(0);
    long long diam = 0;
    int number_of_bad_pairs = 0;
    //cout << t_alpha2 << ' ' << t_omega2 << endl;
    for(int i = 0; i < n2; ++i) {
        if(!is_node[i]) continue;
        vector<long long>t = origin_foremost_path(G2, i, n2, t_alpha2, t_omega2);
        //cout << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << ' ' << t[4] << ' ' << t[5] << endl;
        for(int j = 0; j < n2; ++j) {
            if (is_node[j] and t[j] < oo) {
                if(e[i] < t[j]) {
                    e[i] = t[j];
                    furthest_node[i] = j;
                }
                //e[i] = max(e[i], t[j]);
                //if(i == 1056 and j == 771) cout << t[j] << endl; // OK!
                if (t[j] == oo) number_of_bad_pairs++;
            }
        }

        //cout << e[i] << endl;
    }

    //cout << number_of_bad_pairs << endl;
    //int number_nodes=0;
    //for(int j = 0; j < n; ++j){
    //    if(is_node[j]) number_nodes++;
    //}
    //cout << number_nodes << ' ' << number_of_bad_pairs << endl; // 1547, 1107093 out of 2393209 (46%)


    for(int i = 0; i < n2; ++i) {
        //if(e[i] >= oo) cout << i << endl;
        diam = max(diam, e[i]);
    }
    //cerr << n << endl;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            cerr << par[i][j] << ' ';
        }
        cerr << endl;
    }
    return diam;
}



