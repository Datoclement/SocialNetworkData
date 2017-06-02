#include "relation.hpp"

#include <ctime>
#include <iostream>
#include <string>

#include "mpi.h"

using namespace std;

void performance_test(string dataset,
        string hc_output,string itf_output,string mpi_output,string seq_output){

    int rank,size,root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    vector<relation> rs(3,relation(dataset));
    vector<string> strs;
    strs.push_back("(x,y)");
    strs.push_back("(y,z)");
    strs.push_back("(z,x)");

    if(rank==root) {
        cout << "============================" << endl;
        cout << "\tdataset: " << dataset << endl;
        cout << "\tnumproc: " << size << endl;
    }

    if(rank==root) cout << "hybercube join algorithm..." << endl;
    const clock_t t11 = clock();
    relation hcr = relation::join_hc(rs,strs);
    const clock_t t12 = clock();
    if(rank==root)hcr.save(hc_output);
    const clock_t t13 = clock();
    hcr.free();
    if(rank==root){
        cout << "\ttotl:\t" << (t12-t11)/(double)CLOCKS_PER_SEC << endl;
        cout << "\tsave:\t" << (t13-t12)/(double)CLOCKS_PER_SEC << endl;
    }

    if(rank==root) cout << "instant-transfer join algorithm..." << endl;
    const clock_t t21 = clock();
    relation itfr = relation::join_itf(rs,strs);
    const clock_t t22 = clock();
    if(rank==root)itfr.save(itf_output);
    const clock_t t23 = clock();
    itfr.free();
    if(rank==root){
        cout << "\ttotl:\t" << (t22-t21)/(double)CLOCKS_PER_SEC << endl;
        cout << "\tsave:\t" << (t23-t22)/(double)CLOCKS_PER_SEC << endl;
    }

    if(rank==root) cout << "naive-mpi join algorithm..." << endl;
    const clock_t t31 = clock();
    relation mpir = relation::join_mpi(rs,strs);
    const clock_t t32 = clock();
    if(rank==root)mpir.save(mpi_output);
    const clock_t t33 = clock();
    mpir.free();
    if(rank==root){
        cout << "\ttotl:\t" << (t32-t31)/(double)CLOCKS_PER_SEC << endl;
        cout << "\tsave:\t" << (t33-t32)/(double)CLOCKS_PER_SEC << endl;
    }

    // if(rank==root) cout << "sequential join algorithm..." << endl;
    // const clock_t t41 = clock();
    // relation seqr = relation::join_seq(rs,strs);
    // const clock_t t42 = clock();
    // if(rank==root)seqr.save(mpi_output);
    // const clock_t t43 = clock();
    // seqr.free();
    // if(rank==root){
    //     cout << "\ttotl:\t" << (t42-t41)/(double)CLOCKS_PER_SEC << endl;
    //     cout << "\tsave:\t" << (t43-t42)/(double)CLOCKS_PER_SEC << endl;
    // }

    if(rank==root) {
        cout << "============================" << endl;
    }
}

void debug_test(string hc_output,string itf_output,string mpi_output){
    int rank,root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    vector<relation> rs(3,relation("test"));
    vector<string> strs=*new vector<string>();
    strs.push_back("(x,y)");strs.push_back("(y,z)");strs.push_back("(z,x)");

    relation hcr = relation::join_hc(rs,strs);
    if(rank==root)hcr.save(hc_output);

    relation itfr = relation::join_itf(rs,strs);
    if(rank==root)itfr.save(itf_output);

    relation mpir = relation::join_mpi(rs,strs);
    if(rank==root)mpir.save(mpi_output);

}

int main(int argc,char** argv){

    // string hc_output("res_hc"), itf_output("res_itf"), mpi_output("res_mpi"), seq_output("res_seq");
    string hc_output("res"), itf_output("res"), mpi_output("res"), seq_output("res");

    MPI_Init(&argc,&argv);

    //test for performance
    // string dataset("twitter");
    string dataset("dblp");
    // string dataset("facebook");
    performance_test(dataset,hc_output,itf_output,mpi_output,seq_output);

    // //test for debugging
    // debug_test(hc_output,itf_output,mpi_output,seq_output);

    MPI_Finalize();

    return 0;
}
