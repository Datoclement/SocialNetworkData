#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;

class relation {

private:

    string id;
    int arity,memsize;
    vector<vector<int> > members;
    struct cmp{
        vector<int> permutation;
        cmp(int arity);
        cmp(vector<int> _perm);
        bool operator()(vector<int> a,vector<int> b)const;
    };
    class pattern{
    public:
        int size;
        vector<string> vars;
        map<string,vector<int> > positions;
        pattern(string pat);
        static vector<vector<int> > find_perm(pattern pat1,pattern pat2);
    };
    relation();
    relation(vector<vector<int> > _mem);
    relation filter(pattern pat);

public:

    relation(string str);
    friend ostream& operator<<(ostream& out, const relation& rel);
    void sort();
    void sort(const vector<int>& permutation);
    void random_seed(int seed);
    vector<int> random();
    void save(string filename);

    static relation join(relation r1, relation r2, string pat1, string pat2);
};
