#include "relation.hpp"

#include <iostream>
#include <string>

using namespace std;

int main(int argc,char** argv){

    string id("facebook");
    vector<relation> rs(2,*new relation(id));
    vector<string> ps(2);
    ps[0] = "(x,y)";
    ps[1] = "(y,z)";
    // ps[2] = "(z,x)";
    relation r = relation::join(rs,ps);
    r.save("test.txt");
}
