#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>

#include "relation.hpp"

#define DEBUG 1
#ifdef DEBUG
#define print_array1d(v) {\
    for(int i=0;i<(v).size();i++)\
        cout << (v)[i] << ' ';\
    cout << endl;\
}
#define print_array2d(v) {\
    for(int i=0;i<(v).size();i++){\
        for(int j=0;j<(v)[0].size();j++){\
            cout << (v)[i][j] << ' ';\
        }\
        cout << endl;\
    }\
    cout << endl;\
}
#endif

using namespace std;

relation::relation(){}

relation::relation(vector<vector<int> > mems):
    memsize(mems.size()),
    members(mems){
        if(!mems.empty())arity=mems[0].size();
}

relation::relation(string name):id(name){
    ifstream data("../data/"+name+".txt");
    if(!data.is_open()){
        cout << "Read data failure." << endl;
        return;
    }
    data >> arity >> memsize;
    members = *new vector<vector<int> >(memsize,*new vector<int>(arity));
    for(int i=0;i<memsize;i++)
        for(int j=0;j<arity;j++)
            data >> members[i][j];
}

relation::cmp::cmp(int arity):permutation(arity){
    for(int i=0;i<arity;i++)
        permutation[i] = i;
}

relation::cmp::cmp(vector<int> perm):permutation(perm){}

bool relation::cmp::operator()(vector<int> a,vector<int> b)const{
    int l = a.size();
    for(int i=0;i<l;i++){
        if(a[permutation[i]]<b[permutation[i]])return true;
        if(a[permutation[i]]>b[permutation[i]])return false;
    }
    return false;
}

bool relation::pattern::is_valid_char(char c){
    return (c>='0'&&c<='9')||(c>='a'&&c<='z')||(c>='A'&&c<='Z');
}

relation::pattern::pattern(string pat){
    int l = pat.size();
    if(l==0||pat[0]!='('||pat[l-1]!=')')
        throw "Not valid pattern";
    for(int i=1;i<l-1;i++){
        if(is_valid_char(pat[i]))continue;
        if(pat[i]!=',')
            throw "Not valid pattern";
        if(pat[i]==',')
            if((!is_valid_char(pat[i-1]))||(!is_valid_char(pat[i+1])))
                throw "Not valid pattern";
    }
    int last = 1;
    int cur = 0;
    for(int i=1;i<l;i++){
        if(pat[i]==','||pat[i]==')'){
            string var = pat.substr(last,i-last);
            if(positions[var].size()==0)vars.push_back(var);
            positions[var].push_back(cur);
            cur++;
            last = i+1;
        }
    }
    size = cur;
}

relation::pattern::pattern(vector<string> newvars,map<string,vector<int> > newpositions):
    vars(newvars),positions(newpositions),size(newvars.size()){}

vector<string> relation::pattern::find_comm(pattern pat1,pattern pat2){
    vector<string> communs;
    for(int i=0;i<pat1.vars.size();i++)
        for(int j=0;j<pat2.vars.size();j++)
            if(pat1.vars[i]==pat2.vars[j])
                communs.push_back(pat1.vars[i]);
    return communs;
}

vector<vector<int> > relation::pattern::find_perm(pattern pat1,pattern pat2){
    vector<vector<int> > perms(2);
    perms[0] = *new vector<int>(pat1.vars.size());
    perms[1] = *new vector<int>(pat2.vars.size());

    vector<string> communs = pattern::find_comm(pat1,pat2);
    map<string,bool> is_commun;
    map<string,int> pos_commun;
    for(int i=0;i<communs.size();i++){
        is_commun[communs[i]]=true;
        pos_commun[communs[i]]=i;
    }

    int h,t;
    h=0;t=pat1.vars.size()-1;
    for(int i=0;i<pat1.vars.size();i++){
        if(is_commun[pat1.vars[i]]){
            perms[0][pos_commun[pat1.vars[i]]] = i; h++;
        }
        else {
            perms[0][t] = i; t--;
        }
    }
    h=0;t=pat2.vars.size()-1;
    for(int i=0;i<pat2.vars.size();i++){
        if(is_commun[pat2.vars[i]]){
            perms[1][pos_commun[pat2.vars[i]]] = i; h++;
        }
        else {
            perms[1][t] = i; t--;
        }
    }

    return perms;
}

relation::pattern relation::pattern::join(pattern pat1,pattern pat2){
    vector<string> comm = pattern::find_comm(pat1,pat2);
    vector<vector<int> > perms = pattern::find_perm(pat1,pat2);
    map<string,bool> is_commun;
    for(int i=0;i<comm.size();i++)is_commun[comm[i]]=true;
    vector<string> newvars(pat1.vars.size()+pat2.vars.size()-comm.size());
    map<string,vector<int> > newpositions;
    for(int i=0;i<pat1.vars.size();i++){
        newvars[i] = pat1.vars[i];
        newpositions[(string)pat1.vars[i]].push_back(pat1.positions[(string)pat1.vars[i]][0]);
    }
    vector<string> _positions1rev(perms[1].size());
    for(int i=0;i<pat2.vars.size();i++){
        string curstr = pat2.vars[i];
        for(int j=0;j<pat2.positions[curstr].size();j++){
            int curpos = pat2.positions[curstr][j];
            _positions1rev[curpos] = curstr;
        }
    }
    for(int i=comm.size();i<pat2.vars.size();i++){
        string curstr = _positions1rev[perms[1][i]];
        int curpos = i-comm.size()+pat1.vars.size();
        newvars[curpos] = curstr;
        newpositions[curstr].push_back(curpos);
    }
    return *new pattern(newvars,newpositions);
}

relation relation::filter(pattern pat){
    if(pat.size!=arity)
        throw "Unmatched pattern against relation";
    vector<vector<int> > newmems;
    for(int i=0;i<memsize;i++){
        bool flag = true;
        for(int j=0;j<pat.vars.size();j++){
            vector<int> curpos = pat.positions[(string)pat.vars[j]];
            if(curpos.size()==1)continue;
            for(int k=0;k<curpos.size()-1;k++){
                if(members[i][curpos[k]]!=members[i][curpos[k+1]]){
                    flag = false;
                    break;
                }
            }
            if(!flag)break;
        }
        if(flag){
            vector<int> newmem;
            for(int j=0;j<pat.vars.size();j++){
                newmem.push_back(members[i][pat.positions[(string)pat.vars[j]][0]]);
            }
            newmems.push_back(newmem);
        }
    }
    return *new relation(newmems);
}

void relation::sort(){
    std::sort(members.begin(),members.end(),*new cmp(arity));
}

void relation::sort(const vector<int>& permutation){
    if(permutation.size()!=arity)
        throw "Non valid permutation";
    std::sort(members.begin(),members.end(),*new cmp(permutation));
}

void relation::random_seed(int seed){
    srand(seed);
}

vector<int> relation::random(){
    return members[rand()%members.size()];
}

void relation::save(string filename){
    ofstream file(filename);
    if(file.is_open()){
        string str = "";
        str += to_string(arity) + " " + to_string(members.size()) + "\n";
        for(int i=0;i<members.size();i++){
            for(int j=0;j<arity;j++){
                str += to_string(members[i][j]);
                if(j+1<arity) str += " ";
                else str += "\n";
            }
        }
        file << str;
        file.close();
    }
    else cout << "Saving " << id << " failure." << endl;
}

int cmpf(vector<int> a,vector<int> b,int l,
        const vector<int>& perma,const vector<int>& permb){
    for(int i=0;i<l;i++){
        if(a[perma[i]]<b[permb[i]]){
            return -1;
        }
        if(a[perma[i]]>b[permb[i]]){
            return 1;
        }
    }
    return 0;
}

relation relation::join(relation r1,relation r2,pattern pat1,pattern pat2){
    relation f1 = r1.filter(pat1),f2 = r2.filter(pat2);
    vector<string> comm = pattern::find_comm(pat1,pat2);
    vector<vector<int> > perms = pattern::find_perm(pat1,pat2);
    f1.sort(perms[0]); f2.sort(perms[1]);
    //TODO: the algorithm
    vector<vector<int> > newmems;
    int p1=0,p2=0,l=comm.size();
    cout << "test" << endl;
    while(p1<f1.memsize&&p2<f2.memsize){
        int stat = cmpf(f1.members[p1],f2.members[p2],l,perms[0],perms[1]);
        if(stat > 0) p2++;
        if(stat < 0) p1++;
        if(stat == 0){
            int l1=0,l2=1;
            while(p1+l1<f1.memsize
                    &&cmpf(f1.members[p1],f1.members[p1+l1],l,perms[0],perms[0])==0)l1++;
            while(p2+l2<f2.memsize
                    &&cmpf(f2.members[p2],f2.members[p2+l2],l,perms[1],perms[1])==0)l2++;
            for(int i=0;i<l1;i++){
                for(int j=0;j<l2;j++){
                    vector<int> newmem(f1.members[p1+i]);
                    for(int k=l;k<f2.arity;k++){
                        newmem.push_back(f2.members[p2+j][perms[1][k]]);
                    }
                    newmems.push_back(newmem);
                }
            }
            p1+=l1;
            p2+=l2;
        }
    }
    cout << "test" << endl;
    return *new relation(newmems);
    // return *new relation();
}

relation relation::join(relation r1,relation r2,string patstr1,string patstr2){
    return relation::join(r1,r2,*new pattern(patstr1),*new pattern(patstr2));
}

relation relation::join(vector<relation> rs,vector<string> pats){
    if(rs.size()==1)return rs[0];
    relation res = rs[0];
    pattern pat = *new pattern(pats[0]);
    for(int i=1;i<rs.size();i++){
        res = relation::join(res,rs[i],pat,*new pattern(pats[i]));
        pat = pattern::join(pat,*new pattern(pats[i]));
    }
    return res;
}

// relation relation::multijoin(vector<relation> rs,vector<string> patstrs){
//     if(rs.size()==1)return rs[0];
//     //1.find
// }

ostream& operator<<(ostream &out, relation::pattern &pat){
    out << "#===pattern===#" << endl;
    for(int i=0;i<pat.vars.size();i++){
        out << pat.vars[i] << ": ";
        vector<int> curposs = pat.positions[(pat.vars[i])];
        for(int j=0;j<curposs.size();j++){
            out << curposs[j] << ' ';
        }
        out << endl;
    }
    out << "#=============#" << endl;
    return out;
}

ostream& operator<<(ostream &out, const relation &rel){
    out << "#==begin==#" << endl;
    out << rel.id << " arity:" << rel.arity << endl;
    for(int i=0;i<rel.members.size();i++){
        for(int j=0;j<rel.arity;j++){
            out << rel.members[i][j] << ' ';
        }
        out << endl;
    }
    out << rel.id << " arity:" << rel.arity << endl;
    out << "#==end==#" << endl;
    return out;
}
