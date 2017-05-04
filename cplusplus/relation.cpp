#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>

#include "relation.hpp"

using namespace std;

relation::relation(){}

relation::relation(vector<vector<int> > mems):
    arity(mems[0].size()),
    memsize(mems.size()),
    members(mems){}

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

relation::pattern::pattern(string pat){
    int l = pat.size();
    if(l==0||pat[0]!='('||pat[l-1]!=')')
        throw "Not valid pattern";
    for(int i=1;i<l-1;i++){
        if(pat[i]>='a'&&pat[i]<='z')continue;
        if(pat[i]!=',')
            throw "Not valid pattern";
        if(pat[i]==',')
            if(pat[i-1]<'a'||pat[i-1]>'z'
                ||pat[i+1]<'a'||pat[i+1]>'z')
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

vector<vector<int> > relation::pattern::find_perm(pattern pat1,pattern pat2){
    vector<vector<int> > perms(2);
    perms[0] = *new vector<int>(pat1.size);
    perms[1] = *new vector<int>(pat2.size);

    map<string,bool> is_commun;
    vector<string> communs;
    for(int i=0;i<pat1.vars.size();i++)
        for(int j=0;j<pat2.vars.size();j++)
            if(pat1.vars[i]==pat2.vars[j]){
                communs.push_back(pat1.vars[i]);
                is_commun[pat1.vars[i]]=true;
            }

    int h,t;
    h=0;t=pat1.size-1;
    for(int i=0;i<communs.size();i++){
        vector<int>& cur = pat1.positions[(string)communs[i]];
        perms[0][h] = cur[0]; h++;
        for(int j=1;j<cur.size();j++){
            perms[0][t] = cur[j]; t--;
        }
    }
    for(int i=0;i<pat1.vars.size();i++){
        if(is_commun[pat1.vars[i]])continue;
        vector<int>& cur = pat1.positions[pat1.vars[i]];
        for(int j=0;j<cur.size();j++){
            perms[0][h] = cur[j];h++;
        }
    }
    h=0;t=pat2.size-1;
    for(int i=0;i<communs.size();i++){
        vector<int>& cur = pat2.positions[communs[i]];
        perms[1][h] = cur[0]; h++;
        for(int j=1;j<cur.size();j++){
            perms[1][t] = cur[j]; t--;
        }
    }
    for(int i=0;i<pat2.vars.size();i++){
        if(is_commun[pat2.vars[i]])continue;
        vector<int>& cur = pat2.positions[pat2.vars[i]];
        for(int j=0;j<cur.size();j++){
            perms[1][h] = cur[j];h++;
        }
    }

    return perms;
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
        if(flag) newmems.push_back(members[i]);
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

relation relation::join(relation r1,relation r2,string patstr1,string patstr2){
    pattern pat1(patstr1),pat2(patstr2);
    relation f1 = r1.filter(pat1),f2 = r2.filter(pat2);
    vector<vector<int> > perms = pattern::find_perm(pat1,pat2);
    f1.sort(perms[0]); f2.sort(perms[1]);
    //TODO: the algorithm
    return *new relation();
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
