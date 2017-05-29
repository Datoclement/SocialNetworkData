#include "relation.hpp"

#include <ctime>
#include <iostream>
#include <string>

#include "mpi.h"

using namespace std;

int main(int argc,char** argv){

    const clock_t s = clock();
    int rank,root=0;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(true){
            string source("dblp");
            relation &r1=*new relation(source);
            relation &r2=r1,&r3=r1;
            // if(rank==root){
            //     r1=*new relation(source);
            //     r2=*new relation(source);
        //     r3=*new relation(source);
        // }
        vector<relation> rs;rs.push_back(r1);rs.push_back(r2);rs.push_back(r3);
        vector<string> strs;strs.push_back("(x,y)");strs.push_back("(y,z)");strs.push_back("(z,x)");

        const clock_t t11 = clock();
        if(rank==root) cout << " hc algo..." << endl;
        relation r = relation::join_hypercube(rs,strs);
        const clock_t t12 = clock();
        // r.sort();
        if(rank==root)r.save("res0");
        const clock_t t13 = clock();

        const clock_t t21 = clock();
        if(rank==root) cout << " itf algo..." << endl;
        relation _r = relation::join_mpi(rs,strs);
        const clock_t t22 = clock();
        // _r.sort();
        if(rank==root)_r.save("res1");
        const clock_t t23 = clock();
        if(rank==root){
            cout << (t12-t11)/(double)CLOCKS_PER_SEC << " for hypercube join"<< endl;
            cout << (t13-t12)/(double)CLOCKS_PER_SEC << " for hypercube save"<< endl;
            cout << (t22-t21)/(double)CLOCKS_PER_SEC << " for instant-transfer join"<< endl;
            cout << (t23-t22)/(double)CLOCKS_PER_SEC << " for instant-transfer save"<< endl;

        }
    }
    else{
        string source("test");
        relation r1,r2,r3;
        if(rank==root){
            r1=*new relation(source);
            r2=*new relation(source);
            r3=*new relation(source);
        }
        vector<relation> rs=*new vector<relation>();
        rs.push_back(r1);rs.push_back(r2);rs.push_back(r3);
        vector<string> strs=*new vector<string>();
        strs.push_back("(x,y)");strs.push_back("(y,z)");strs.push_back("(z,x)");
        // relation r = relation::join_mpi(r1,r2,"(x,y,z)","(y,x,h)");
        relation r = relation::join_hypercube(rs,strs);
        if(rank==root)r.save("res3");
        relation _r = relation::join_mpi(rs,strs);
        if(rank==root)r.save("res4");
    }
    const clock_t t3 = clock();

    MPI_Finalize();
    const clock_t e = clock();


}
