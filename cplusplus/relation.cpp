#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>

#include "relation.hpp"

using namespace std;

relation::relation(string name):arity(2),id(name){
    while(true){
        vector<int> newmem(arity,0);
        bool flag = false;
        for(int i=0;i<arity;i++)
            if(!(cin >> newmem[i])){
                flag = true;
                break;
            }
        if (flag) break;
        members.push_back(newmem);
    }
}

relation::cmp::cmp(int arity):permutation(arity){
    for(int i=0;i<arity;i++)
        permutation[i] = i;
}

relation::cmp::cmp(vector<int> _perm):permutation(_perm){}

bool relation::cmp::operator()(vector<int> a,vector<int> b)const{
    int l = a.size();
    for(int i=0;i<l;i++){
        if(a[permutation[i]]<b[permutation[i]])return true;
        if(a[permutation[i]]>b[permutation[i]])return false;
    }
    return false;
}

void relation::sort(){
    std::sort(members.begin(),members.end(),*new cmp(arity));
}

void relation::sort(const vector<int> permutation){
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
