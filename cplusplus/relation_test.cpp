#include "relation.hpp"

#include <iostream>
#include <string>

using namespace std;

int main(int argc,char** argv){

    // string id(argv[0]);
    relation r("facebook");
    vector<int> perm(2,0);
    perm[0]=1;
    r.sort(perm);
    cout << r << endl;
    r.random_seed(0);
    for(int i=0;i<10;i++){
        vector<int> cur = r.random();
        for(int j=0;j<2;j++){
            cout << cur[j] << ' ';
        }
        cout << endl;
    }
    r.save("facebook_sorted_by_second_attribute");
}
