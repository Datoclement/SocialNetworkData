#include <iostream>
#include <vector>
#include <string>

using namespace std;

class relation {

private:
    string id;
    int arity;
    vector<vector<int> > members;
    struct cmp{
        vector<int> permutation;
        cmp(int arity);
        cmp(vector<int> _perm);
        bool operator()(vector<int> a,vector<int> b)const;
    };
public:
    relation(string str);
    friend ostream& operator<<(ostream& out, const relation& rel);
    void sort();
    void sort(const vector<int> permutation);
    void random_seed(int seed);
    vector<int> random();
    void save(string filename);
};
