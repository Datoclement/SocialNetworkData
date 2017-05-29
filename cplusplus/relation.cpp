#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <iterator>
#include <utility>
#include <ctime>

#include "mpi.h"

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

// #define POST_TRANSFER true

relation::relation():arity(0),memsize(0){}

relation::relation(vector<vector<int> >& mems):
    memsize(mems.size()),
    members(mems){
        if(!mems.empty())arity=mems[0].size();
        else arity=0;
}

relation::relation(string name):id(name){
    ifstream data(("../data/"+name+".txt").c_str());
    // ifstream data;
    // data.open("../data/"+name+".txt",ifstream::in);
    if(!data.is_open()){
        cout << "Read data failure." << endl;
        return;
    }
    data >> arity >> memsize;

    members = vector<vector<int> >(memsize,vector<int>(arity));
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

vector<string> relation::pattern::find_comm(vector<pattern> pats){
    vector<string> res;
    map<string,bool> exist;
    map<string,bool> inres;
    for(int i=0;i<pats.size();i++){
        for(int j=0;j<pats[i].vars.size();j++){
            string& cur=pats[i].vars[j];
            if(!exist[cur]){
                exist[cur]=true;continue;
            }
            if(inres[cur])continue;
            res.push_back(cur);
            inres[cur]=true;
        }
    }
    return res;
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
    map<string,bool> is_commun;
    for(int i=0;i<comm.size();i++)is_commun[comm[i]]=true;
    vector<string> newvars(pat1.vars);
    map<string,vector<int> > newpositions;
    for(int i=0;i<pat1.vars.size();i++)newpositions[pat1.vars[i]].push_back(i);
    int cnt = pat1.vars.size();
    for(int i=0;i<pat2.vars.size();i++){
        string curstr = pat2.vars[i];
        if(is_commun[curstr])continue;
        newvars.push_back(pat2.vars[i]);
        newpositions[pat2.vars[i]].push_back(cnt);
        cnt++;
    }
    return *new pattern(newvars,newpositions);
}

relation::cubehash::cubehash(vector<int>_ms,vector<string>&comm,relation::pattern pat):ms(_ms){

    //get number of processors
    int ttp=1;
    for(int i=0;i<ms.size();i++)ttp*=ms[i];

    //get is_patvars function(structure)
    map<string,bool> _is_patvars;
    for(int i=0;i<pat.vars.size();i++)_is_patvars[pat.vars[i]]=true;

    //get the position if the var exists in pat
    varspos=*new vector<int>(ms.size(),-1);
    for(int i=0;i<ms.size();i++){
        if(_is_patvars[comm[i]]){
            varspos[i]=pat.positions[comm[i]][0];
        }
    }

    //get number of possible hash values
    int vlen=1;
    for(int i=0;i<ms.size();i++){
        if(_is_patvars[comm[i]])vlen*=ms[i];
    }

    //put each processor in the list of correspond hash value
    values=*new vector<vector<int> >(vlen);
    for(int i=0;i<ttp;i++){
        int curl=ttp,cnt=0,cur=i;
        for(int j=0;j<ms.size();j++){
            curl /= ms[j];
            if(_is_patvars[comm[j]]){
                cnt*=ms[j];
                cnt+=cur/curl;
            }
            cur %= curl;
        }
        values[cnt].push_back(i);
    }
}

vector<int>& relation::cubehash::get_value(vector<int>& v){
    int id=0;
    for(int i=0;i<ms.size();i++){
        if(varspos[i]>=0){
            id*=ms[i];
            id+=v[varspos[i]]%ms[i];
        }
    }
    return values[id];
}

pair<relation,relation::pattern> relation::filter(pattern pat){
    map<string,vector<int> > newpositions;
    for(int i=0;i<pat.vars.size();i++){
        newpositions[pat.vars[i]].push_back(i);
    }
    sort();
    vector<vector<int> > newmems;
    for(int i=0;i<memsize;i++){
        if(i>0&&members[i]==members[i-1])continue;
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

    return pair<relation,relation::pattern>(* new relation(newmems),pattern(pat.vars,newpositions));
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
    cout << "saving to " << filename << "..." << endl;
    ofstream file;
    file.open(("../output/"+filename).c_str());
    if(file.is_open()){
        file << arity << " " << members.size() << "\n";
        for(int i=0;i<members.size();i++){
            for(int j=0;j<arity;j++){
                file << (members[i][j]);
                if(j+1<arity) file << " ";
                else file << "\n";
            }
        }
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

//sequential join
relation& relation::join(relation& r1,relation& r2,pattern _pat1,pattern _pat2){
    if(r1.memsize==0||r2.memsize==0)return * new relation();
    pair<relation,relation::pattern> pa1 = r1.filter(_pat1),pa2 = r2.filter(_pat2);
    relation& f1 = pa1.first, f2 = pa2.first;
    pattern pat1 = pa1.second, pat2 = pa2.second;
    vector<string> comm = pattern::find_comm(pat1,pat2);
    vector<vector<int> > perms = pattern::find_perm(pat1,pat2);
    map<string,bool> is_commun;
    for(int i=0;i<comm.size();i++)is_commun[comm[i]]=true;

    pattern p = pattern::join(pat1,pat2);
    f1.sort(perms[0]); f2.sort(perms[1]);
    //TODO: the algorithm
    vector<vector<int> > newmems;
    int p1=0,p2=0,l=comm.size();
    while(p1<f1.memsize&&p2<f2.memsize){
        int stat = cmpf(f1.members[p1],f2.members[p2],l,perms[0],perms[1]);
        if(stat > 0) p2++;
        if(stat < 0) p1++;
        if(stat == 0){
            int l1=0,l2=0;
            while(p1+l1<f1.memsize
                    &&cmpf(f1.members[p1],f1.members[p1+l1],l,perms[0],perms[0])==0)l1++;
            while(p2+l2<f2.memsize
                    &&cmpf(f2.members[p2],f2.members[p2+l2],l,perms[1],perms[1])==0)l2++;
            for(int i=0;i<l1;i++){
                for(int j=0;j<l2;j++){
                    vector<int> newmem(f1.members[p1+i]);
                    for(int k=0;k<pat2.vars.size();k++){
                        if(is_commun[pat2.vars[k]])continue;
                        newmem.push_back(f2.members[p2+j][k]);
                    }
                    newmems.push_back(newmem);
                }
            }
            p1+=l1;
            p2+=l2;
        }
    }
    return *new relation(newmems);
}

relation& relation::join(relation& r1,relation& r2,string patstr1,string patstr2){
    return relation::join(r1,r2,*new pattern(patstr1),*new pattern(patstr2));
}

int _hash(int a,int b){
    return a%b;
}

relation& relation::join_mpi(relation& r1,relation& r2,string patstr1,string patstr2){
    return relation::join_mpi_hash(r1,r2,*new pattern(patstr1),*new pattern(patstr2),_hash);

}

//naive parallel join
relation& relation::join_mpi_hash(relation& r1,relation& r2,
        pattern _pat1, pattern _pat2, hashtype* hash){
    int id, size; int root = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    pair<relation,relation::pattern> p1 = r1.filter(_pat1),p2 = r2.filter(_pat2);
    relation &f1 = p1.first, &f2 = p2.first;
    pattern pat1 = p1.second, pat2 = p2.second;
    vector<string> comm = pattern::find_comm(pat1,pat2);
    vector<vector<int> > mem1;
    vector<vector<int> > mem2;
    int pos1=pat1.positions[comm[0]][0];
    int pos2=pat2.positions[comm[0]][0];

    for (int i=0; i<f1.memsize;i++){
        if (hash(f1.members[i][pos1],size)==id){
            mem1.push_back(f1.members[i]) ;
        }
    }
    for (int i=0; i<f2.memsize;i++){
        if (hash(f2.members[i][pos2],size)==id){
            mem2.push_back(f2.members[i]);
        }
    }
    relation &tmp1=*new relation(mem1),&tmp2=*new relation(mem2);
    relation &local= join(tmp1,tmp2,pat1,pat2);

    if(id!=root){
        MPI_Send(&local.memsize,1,MPI_INT,root,0,MPI_COMM_WORLD);
        if(local.memsize>0){
            MPI_Send(&local.arity,1,MPI_INT,root,1,MPI_COMM_WORLD);
            for(int i=0;i<local.memsize;i++){
                int* toSend=&(local.members[i][0]);
                MPI_Send(toSend,local.members[i].size(),MPI_INT,root,i+2,MPI_COMM_WORLD);
            }
        }
        return *new relation();
    }
    else{
        vector<vector<int> > mem;
        mem.insert(mem.end(),local.members.begin(),local.members.end());
        for(int i=1;i<size;i++){
            int len;
            MPI_Recv(&len,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(len>0){
                int ar;
                MPI_Recv(&ar,1,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                int* tmp = new int[ar];
                for(int j=0;j<len;j++){
                    MPI_Recv(tmp,ar,MPI_INT,i,j+2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    mem.push_back(vector<int>(tmp,tmp+ar));
                }
                delete [] tmp;
            }
        }

        return *new relation(mem);
    }
}

relation& relation::join(vector<relation>& rs,vector<pattern>& pats){
    if(rs.size()==0)return *new relation();
    relation *res=&rs[0];
    pattern pat = pats[0];
    for(int i=1;i<rs.size();i++){
        *res = join(*res,rs[i],pat,pats[i]);
        pat = pattern::join(pat,pats[i]);
    }
    return *res;
}

//naive multijoin
relation& relation::join(vector<relation>& rs,vector<string>& patstrs){
    if(rs.size()==1)return *new relation();
    vector<pattern> pats;
    for(int i=0;i<patstrs.size();i++){
        pats.push_back(pattern(patstrs[i]));
    }
    return join(rs,pats);
}

//distribute members according to a single variable
vector<vector<int> >& relation::distribute(vector<vector<int> >& mems,int pos){
    int rank,size;int root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //get ar
    int ar=0;if(mems.size()>0)ar=mems[0].size();
    int* arv=new int[size];
    MPI_Gather(&ar,1,MPI_INT,arv,1,MPI_INT,root,MPI_COMM_WORLD);
    if(rank==root){
        for(int i=0;i<size;i++){
            if(arv[i]>0){
                ar=arv[i];
                break;
            }
        }
    }
    MPI_Bcast(&ar,1,MPI_INT,root,MPI_COMM_WORLD);
    delete [] arv;
    if(ar==0)return *new vector<vector<int> >();
    //get send lens and send displ
    int* slen=new int[size];
    for(int i=0;i<size;i++)slen[i]=0;
    for(int i=0;i<mems.size();i++)slen[_hash(mems[i][pos],size)]+=ar;
    int* sdis=new int[size];
    sdis[0]=0;
    for(int i=1;i<size;i++)sdis[i]=sdis[i-1]+slen[i-1];

    //get recv lens and recv displ
    int* rlen=new int[size];
    int* rdis=new int[size];
    MPI_Alltoall(slen,1,MPI_INT,rlen,1,MPI_INT,MPI_COMM_WORLD);
    rdis[0]=0;
    for(int i=1;i<size;i++)rdis[i]=rdis[i-1]+rlen[i-1];

    //slim data
    int* ipos=new int[size];
    for(int i=0;i<size;i++)ipos[i]=0;
    int* data=new int[mems.size()*ar];
    for(int i=0;i<mems.size();i++){
        int curid=_hash(mems[i][pos],size);
        int curpos=ipos[curid]+sdis[curid];
        ipos[curid]+=ar;
        for(int j=0;j<ar;j++){
            data[curpos+j]=mems[i][j];
        }
    }
    delete [] ipos;

    int ttl=rlen[size-1]+rdis[size-1];
    int* collect=new int[ttl];
    MPI_Alltoallv(data,slen,sdis,MPI_INT,collect,rlen,rdis,MPI_INT,MPI_COMM_WORLD);

    delete [] data;
    delete [] rlen;
    delete [] slen;
    delete [] rdis;
    delete [] sdis;

    vector<vector<int> >& newmems=*new vector<vector<int> >(ttl/ar,vector<int>(ar));
    for(int i=0;i<ttl/ar;i++){
        for(int j=0;j<ar;j++){
            newmems[i][j]=collect[i*ar+j];
        }
    }
    return newmems;

}

//distribute members according to the commun variables tuple
vector<vector<int> >& relation::distribute_cube(vector<vector<int> >& _mems,int source,cubehash cbh){

    int rank,size;int root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);


    int l;if(rank==source)l=_mems.size();
    MPI_Bcast(&l,1,MPI_INT,source,MPI_COMM_WORLD);
    int ar;if(rank==source)ar=_mems[0].size();
    MPI_Bcast(&ar,1,MPI_INT,source,MPI_COMM_WORLD);


    if(l<=2000000){
        int* lenv=new int[size];
        for(int i=0;i<size;i++)lenv[i]=0;
        int* thresv=new int[size];
        // if(rank==source){
        //     for(int i=0;i<size;i++){
        //         cout << lenv[i] <<' ';
        //     }
        //     cout << endl;
        // }

        if(rank==source){
            for(int i=0;i<_mems.size();i++){
                vector<int> cur=cbh.get_value(_mems[i]);
                // print_array1d(cur);
                for(int j=0;j<cur.size();j++){
                    lenv[cur[j]]+=ar;
                }
            }
            thresv[0]=0;
            for(int i=1;i<size;i++)thresv[i]=thresv[i-1]+lenv[i-1];
        }
        // cout << ar << endl;
        // if(rank==source){
        //     for(int i=0;i<size;i++){
        //         cout << lenv[i] <<' ';
        //     }
        //     cout << endl;
        // }
        int len;
        MPI_Scatter(lenv,1,MPI_INT,&len,1,MPI_INT,source,MPI_COMM_WORLD);

        int* pos=new int[size];
        for(int i=0;i<size;i++)pos[i]=0;


        int* data;
        if(rank==source){
            int ttl=thresv[size-1]+lenv[size-1];
            data=new int[ttl];
            for(int i=0;i<_mems.size();i++){
                vector<int> cur=cbh.get_value(_mems[i]);
                for(int j=0;j<cur.size();j++){
                    int curid=cur[j];
                    int curpos=thresv[curid]+pos[curid];
                    pos[curid]+=ar;
                    for(int k=0;k<ar;k++){
                        data[curpos+k]=_mems[i][k];
                    }
                }
            }
        }

        int* collect=new int[len];

        MPI_Scatterv(data,lenv,thresv,MPI_INT,collect,len,MPI_INT,source,MPI_COMM_WORLD);
        vector<vector<int> >& newmems=*new vector<vector<int> >(len/ar,vector<int>(ar));
        for(int i=0;i<len/ar;i++){
            for(int j=0;j<ar;j++){
                newmems[i][j]=collect[i*ar+j];
            }
        }
        // cout << "test" << endl;

        return newmems;
    }
    else{
        int b=(rank<l%size);
        int thres=b*rank*(l/size+1)+(1-b)*((l%size)*(l/size+1)+(rank-l%size)*(l/size));
        int len=l/size+b;
        thres*=ar;len*=ar;

        int* thresv=new int[size];
        int* lenv=new int[size];
        MPI_Gather(&thres,1,MPI_INT,thresv,1,MPI_INT,source,MPI_COMM_WORLD);
        MPI_Gather(&len,1,MPI_INT,lenv,1,MPI_INT,source,MPI_COMM_WORLD);

        int* _memsv=new int[l*ar];
        if(rank==source){
            for(int i=0;i<l;i++){
                for(int j=0;j<ar;j++){
                    _memsv[i*ar+j]=_mems[i][j];
                }
            }
        }

        int* mems_loc=new int[len];
        MPI_Scatterv(_memsv,lenv,thresv,MPI_INT,mems_loc,len,MPI_INT,source,MPI_COMM_WORLD);

        delete [] thresv;
        delete [] lenv;
        delete [] _memsv;

        vector<vector<int> >mems(len/ar,vector<int>(ar));
        for(int i=0;i<len/ar;i++){
            for(int j=0;j<ar;j++){
                mems[i][j]=mems_loc[i*ar+j];
            }
        }
        delete [] mems_loc;
        // print_array2d(mems);
        //get send lens and send displ
        int* slen=new int[size];
        for(int i=0;i<size;i++)slen[i]=0;
        for(int i=0;i<mems.size();i++){
            vector<int> cur=cbh.get_value(mems[i]);
            for(int j=0;j<cur.size();j++){
                slen[cur[j]]+=ar;
            }
        }
        int* sdis=new int[size];
        sdis[0]=0;
        for(int i=1;i<size;i++)sdis[i]=sdis[i-1]+slen[i-1];

        //get recv lens and recv displ
        int* rlen=new int[size];
        int* rdis=new int[size];
        MPI_Alltoall(slen,1,MPI_INT,rlen,1,MPI_INT,MPI_COMM_WORLD);
        rdis[0]=0;
        for(int i=1;i<size;i++)rdis[i]=rdis[i-1]+rlen[i-1];

        //slim data
        int* ipos=new int[size];
        for(int i=0;i<size;i++)ipos[i]=0;
        int* data=new int[slen[size-1]+sdis[size-1]];
        for(int i=0;i<mems.size();i++){
            vector<int> cur=cbh.get_value(mems[i]);
            for(int j=0;j<cur.size();j++){
                int curid=cur[j];
                int curpos=ipos[curid]+sdis[curid];
                ipos[curid]+=ar;
                for(int k=0;k<ar;k++){
                    data[curpos+k]=mems[i][k];
                }
            }
        }
        delete [] ipos;

        // cout << "_f" << endl;
        // cout << "data" << endl;
        // for(int i=0;i<mems.size();i++){
        //     for(int j=0;j<ar;j++){
        //         cout << data[i*ar+j] << ' ';
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        int ttl=rlen[size-1]+rdis[size-1];
        // cout << ttl << endl;
        // for(int i=0;i<size;i++){
        //     cout << "(" << rlen[i] << "," << rdis[i] << "," << slen[i] << "," << sdis[i] << ")" << endl;
        // }
        int* collect=new int[ttl];

        MPI_Alltoallv(data,slen,sdis,MPI_INT,collect,rlen,rdis,MPI_INT,MPI_COMM_WORLD);
        // cout << "test" << endl;


        delete [] data;
        delete [] rlen;
        delete [] slen;
        delete [] rdis;
        delete [] sdis;
        // cout << ttl << "," << ar << endl;
        vector<vector<int> >& newmems=*new vector<vector<int> >(ttl/ar,vector<int>(ar));
        for(int i=0;i<ttl/ar;i++){
            for(int j=0;j<ar;j++){
                newmems[i][j]=collect[i*ar+j];
            }
        }
        delete [] collect;
        return newmems;
    }
}

//merge members into a single relation
relation& relation::merge_vector_to_relation(vector<vector<int> >&mems,int root){
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // cout << rank << " merging..." << endl;
    int len=mems.size();
    int ar=0;if(mems.size()>0)ar=mems[0].size();
    int* lenv=new int[size];
    int* arv=new int[size];
    MPI_Gather(&len,1,MPI_INT,lenv,1,MPI_INT,root,MPI_COMM_WORLD);
    MPI_Gather(&ar,1,MPI_INT,arv,1,MPI_INT,root,MPI_COMM_WORLD);

    if(rank==root){
        for(int i=0;i<size;i++){
            if(arv[i]>0){
                ar=arv[i];
                break;
            }
        }
    }
    MPI_Bcast(&ar,1,MPI_INT,root,MPI_COMM_WORLD);

    int ttl=0;
    if(rank==root){
        for(int i=0;i<size;i++)ttl+=lenv[i];
    }

    int* data=new int[1];
    if(len>0)data=new int[len*ar];
    for(int i=0;i<len;i++){
        for(int j=0;j<ar;j++){
            data[i*ar+j]=mems[i][j];
        }
    }
    int* displ=new int[size];
    if(rank==root){
        displ[0]=0;
        for(int i=1;i<size;i++)displ[i]=displ[i-1]+lenv[i-1]*ar;
    }
    int* res_vect=NULL;
    if(rank==root){
        res_vect=new int[ttl*ar];
    }
    if(rank==root){
        for(int i=0;i<size;i++)lenv[i]*=ar;
    }
    // cout << rank << " waits for gather" << endl;

    MPI_Gatherv(data,len*ar,MPI_INT,res_vect,lenv,displ,MPI_INT,root,MPI_COMM_WORLD);
    delete [] lenv;
    delete [] arv;
    delete [] data;
    delete [] displ;

    // cout << rank << " through gather" << endl;

    if(rank!=root) {
        return *new relation();
    }
    else{
        // cout << "root arrives here" << endl;
        vector<vector<int> > newmems(ttl,vector<int>(ar));
        for(int i=0;i<ttl;i++){
            for(int j=0;j<ar;j++){
                newmems[i][j]=res_vect[i*ar+j];
            }
        }

        // cout << ttl << "," << ar << endl;
        // cout << ttl*ar << endl;
        // vector<vector<int> > newmems=*new vector<vector<int> >(ttl,*new vector<int>(ar));
        // for(int i=0;i<ttl;i++){
        //     cout << "working on " << i << " element" << endl;
        //     for(int j=0;j<ar;j++){
        //         newmems[i][j]=res_vect[i*ar+j];
        //         // newmems[i].push_back(res_vect[i*ar+j]);
        //     }
        // }
        // ofstream file;
        // string str="../data/tmp.txt";
        // file.open(str.c_str());
        // if(file.is_open()){
        //     file << ar << " " << ttl << "\n";
        //     for(int i=0;i<ttl;i++){
        //         for(int j=0;j<ar;j++){
        //             file << res_vect[i*ar+j];
        //             if(j+1==ar) file << "\n";
        //             else file << " ";
        //         }
        //     }
        //     file.close();
        //     cout << "Save tmp complete." << endl;
        // }
        // else cout << "Saving tmp.txt failure." << endl;
        delete [] res_vect;
        return *new relation(newmems);
        // return *new relation("tmp");
    }
}

//instant transfer join
relation& relation::join_mpi(vector<relation>& rs,vector<pattern>& pats,int root){

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    bool quick_return=false;
    if(rank==root&&rs.size()==0){
        quick_return=true;
        MPI_Bcast(&quick_return,1,MPI_INT,root,MPI_COMM_WORLD);
    }
    if(quick_return)return *new relation();
    // if(rank==root)cout << rank << " instant transfer joining..." << endl;
    pair<relation,relation::pattern> p=rs[0].filter(pats[0]);
    vector<vector<int> >* currel=(rank==root?&p.first.members:new vector<vector<int> >());
    pattern curpat=p.second;

    double dt=0,jt=0,mt=0;

    for(int i=1;i<rs.size();i++){
        pair<relation,relation::pattern> p=rs[i].filter(pats[i]);
        vector<string> comm=pattern::find_comm(curpat,p.second);

        int pos1=curpat.positions[comm[0]][0];
        int pos2=p.second.positions[comm[0]][0];

        vector<vector<int> >&nextrel=(rank==root?p.first.members:*new vector<vector<int> >());

        const clock_t s1=clock();
        vector<vector<int> >&currel_local=distribute(*currel,pos1);
        vector<vector<int> >&next_local=distribute(nextrel,pos2);
        // vector<vector<int> >&currel_local=*new vector<vector<int> >();
        // vector<vector<int> >&next_local=*new vector<vector<int> >();
        const clock_t t1=clock();
        dt+=(t1-s1)/(double)CLOCKS_PER_SEC;

        //
        //
        // for(int j=0;j<size;j++){
        //     vector<vector<int> >& tmp1=*new vector<vector<int> >();
        //     tmp1=distribute(*currel,j,pos1);
        //     currel_local.insert(currel_local.end(),tmp1.begin(),tmp1.end());
        //     vector<vector<int> >& tmp2=*new vector<vector<int> >();
        //     tmp2=distribute(nextrel,j,pos2);
        //     next_local.insert(next_local.end(),tmp2.begin(),tmp2.end());
        //
        // }
        const clock_t s2=clock();
        *currel=join(*new relation(currel_local),*new relation(next_local),curpat,pats[i]).members;
        const clock_t t2=clock();
        jt+=(t2-s2)/(double)CLOCKS_PER_SEC;
        // cout << rank << " get result size " << currel->size()
                // << " for " << root << " in " << i << "-th round..." << endl;
        curpat=relation::pattern::join(curpat,p.second);
    }
    const clock_t s3=clock();
    // cout << rank << " entering the merge function..." << endl;
    relation& res=merge_vector_to_relation(*currel,root);
    // if(rank==root)cout << root << " get partial result " << res << endl;
    // if(rank==root)cout << rank << " instant transfer join finished." << endl;
    const clock_t t3=clock();
    mt+=(t3-s3)/(double)CLOCKS_PER_SEC;
    if(rank==root)cout << " itf distribute time: " << dt << endl;
    if(rank==root)cout << " itf join time: " << jt << endl;
    if(rank==root)cout << " itf merge time: " << mt << endl;
    return res;
}

relation& relation::join_mpi(vector<relation>&rs,vector<string>& strs){
    vector<pattern> pats;
    for(int i=0;i<strs.size();i++)
        pats.push_back(pattern(strs[i]));
    return relation::join_mpi(rs,pats,0);
}

vector<int> get_distribution(int size,int vars){
    vector<int> res(vars);
    for(int i=0;i<vars-1;i++){
        res[i]=1;
    }
    res[vars-1]=size;
    bool change=true;
    while(true){
        if(!change)break;
        change=false;
        for(int i=vars-1;i>0;i--){
            for(int j=i-1;j>=0;j--){
                int mxd=res[j];
                for(int k=res[j]+1;k*k<=res[i]*res[j];k++){
                    if((res[i]*res[j])%k==0){
                        mxd=k;
                    }
                }
                if(mxd!=res[j]){
                    change=true;
                    int tmp=res[i]*res[j];
                    res[j]=mxd;
                    res[i]=tmp/mxd;
                }
            }
        }
    }
    return res;
}

relation& relation::join_hypercube(vector<relation>&rs,vector<pattern>& pats,int root){
    if(rs.size()==0)return *new relation();

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    // cout << rank << " begin..." << endl;

    vector<string> comms=pattern::find_comm(pats);
    vector<int> ms=get_distribution(size,comms.size());

    // if(rank==root){
    //     // cout << "variable order: " << endl;
    //     print_array1d(ms);
    //     print_array1d(comms);
    // }
    vector<cubehash> cbh;
    for(int i=0;i<pats.size();i++){
        cbh.push_back(*new cubehash(ms,comms,pats[i]));
    }
    vector<relation>& rs_loc=*new vector<relation>();

    const clock_t s1=clock();
    for(int i=0;i<rs.size();i++){
        // cout << rank << " hypercube-distributing " << i << " ..." << endl;
        // cout << rank << " " << rs[0].memsize << " ..." << endl;
        rs_loc.push_back(*new relation(distribute_cube(rs[i].members,root,cbh[i])));
    }
    const clock_t t1=clock();
    if(rank==root)cout << " hc distribute time: " << (t1-s1)/(double)CLOCKS_PER_SEC << endl;
    // cout << rank << " finish hypercube-distribution." << endl;
    vector<vector<int> > res_loc;

    //hypeflag
    const clock_t s2=clock();
    #ifdef POST_TRANSFER
        // using instant tranfer for each step
        // cout << rank << " instant-transfer-joining..." << endl;
        for(int i=0;i<size;i++){
            relation& tmp=join_mpi(rs_loc,pats,i);
            res_loc.insert(res_loc.end(),tmp.members.begin(),tmp.members.end());
        }
    #else
        //naive join
        // cout << rank << " local joining..." << endl;
        res_loc=(join(rs_loc,pats)).members;
    #endif
    const clock_t t2=clock();
    if(rank==root)cout << " hc join time: " << (t2-s2)/(double)CLOCKS_PER_SEC << endl;

    // cout << rank << " get result size " << res_loc.size() << endl;

    // cout << rank << " enters the merging function..." << endl;
    const clock_t s3=clock();
    relation& res=merge_vector_to_relation(res_loc,root);
    const clock_t t3=clock();
    if(rank==root)cout << " hc merge time: " << (t3-s3)/(double)CLOCKS_PER_SEC << endl;

    return res;

}

relation& relation::join_hypercube(vector<relation>&rs,vector<string>& patstrs){
    vector<pattern> pats;
    for(int i=0;i<patstrs.size();i++){
        pats.push_back(pattern(patstrs[i]));
    }
    return join_hypercube(rs,pats,0);
}

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
