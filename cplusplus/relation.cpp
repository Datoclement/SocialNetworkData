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

// #define FLTR

using namespace std;

//constructor
relation::relation():arity(0),memsize(0),members(NULL){}
relation::relation(int* mems,int a,int m):arity(a),memsize(m),members(mems){}
relation::relation(string name){
    ifstream data(("../data/"+name+".txt").c_str());
    if(!data.is_open()){
        cout << "Read data failure." << endl;
        return;
    }
    data >> arity >> memsize;
    members=new int[arity*memsize];
    for(int i=0;i<memsize;i++)
        for(int j=0;j<arity;j++)
            data >> members[i*arity+j];
}

//cmp class
relation::cmp::cmp(int arity):permutation(arity){
    for(int i=0;i<arity;i++)
        permutation[i] = i;
}
relation::cmp::cmp(vector<int> perm):permutation(perm){}
bool relation::cmp::operator()(int* a,int* b)const{
    int l = permutation.size();
    for(int i=0;i<l;i++){
        if(a[permutation[i]]<b[permutation[i]])return true;
        if(a[permutation[i]]>b[permutation[i]])return false;
    }
    return false;
}

//pattern class
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

//hash_hc class
relation::hash_hc::hash_hc(vector<int>_ms,vector<string>&comm,relation::pattern pat):ms(_ms){

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
vector<int>& relation::hash_hc::get_value(int* v){
    int id=0;
    for(int i=0;i<ms.size();i++){
        if(varspos[i]>=0){
            id*=ms[i];
            id+=v[varspos[i]]%ms[i];
        }
    }
    return values[id];
}

//hash_itf class
relation::hash_itf::hash_itf(int _sz):sz(_sz),r(_sz),pos(0){}
void relation::hash_itf::set_pos(int _pos){pos=_pos;}
void relation::hash_itf::evolve(){r++;}
int relation::hash_itf::get_value(int* key){return key[pos]%r%sz;}

//clean data function
pair<relation,relation::pattern> relation::filter(pattern pat){
    map<string,vector<int> > newpositions;
    for(int i=0;i<pat.vars.size();i++){
        newpositions[pat.vars[i]].push_back(i);
    }
    sort();
    vector<int> idset;
    for(int i=0;i<memsize;i++){
        if(i>0){
            bool equal=true;
            for(int j=0;j<arity;j++){
                if(get(i)[j]!=get(i-1)[j])equal=false;
            }
            if(equal)continue;
        }
        bool flag = true;
        int* cur=get(i);
        for(int j=0;j<pat.vars.size();j++){
            vector<int> curpos = pat.positions[(string)pat.vars[j]];
            if(curpos.size()==1)continue;

            for(int k=0;k<curpos.size()-1;k++){
                if(cur[curpos[k]]!=cur[curpos[k+1]]){
                    flag = false;
                    break;
                }
            }
            if(!flag)break;
        }
        if(flag)idset.push_back(i);
    }
    int newsize=idset.size();
    int* newmems=new int[newsize*arity];
    for(int i=0;i<newsize;i++){
        for(int j=0;j<arity;j++){
            newmems[i*arity+j]=members[idset[i]*arity+j];
        }
    }

    return pair<relation,relation::pattern>(
            * new relation(newmems,arity,newsize),pattern(pat.vars,newpositions));
}

//sort functions
void relation::sort(cmp f){
    vector<int*> heads(memsize);
    for(int i=0;i<memsize;i++){
        heads[i]=members+i*arity;
    }
    std::sort(heads.begin(),heads.end(),f);
    int* newmems=new int[arity*memsize];
    for(int i=0;i<memsize;i++){
        for(int j=0;j<arity;j++){
            newmems[i*arity+j]=heads[i][j];
        }
    }
    //TODO
    members = newmems;
}
void relation::sort(){
    sort(cmp(arity));
}
void relation::sort(const vector<int>& permutation){
    sort(cmp(permutation));
}

//access functions
int* relation::get(int id){return members+id*arity;}
void relation::random_seed(int seed){
    srand(seed);
}
int* relation::random(){
    return members+(rand()%memsize*arity);
}

//helpers (important)
int relation::cmpf(int* a,int* b,int l,vector<int>& perma,vector<int>& permb){
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
relation& relation::merge(relation& r,int root){
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // cout << rank << " merging..." << endl;
    int sz=r.memsize;
    int ar=r.arity;
    int* szv=new int[size];
    int* arv=new int[size];
    MPI_Gather(&sz,1,MPI_INT,szv,1,MPI_INT,root,MPI_COMM_WORLD);
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

    if(ar==0)return *new relation();

    int nsz=0;
    if(rank==root){
        for(int i=0;i<size;i++)nsz+=szv[i];
    }

    int* displ=new int[size];
    if(rank==root){
        displ[0]=0;
        for(int i=1;i<size;i++)displ[i]=displ[i-1]+szv[i-1]*ar;
    }
    int* res_vect=NULL;
    if(rank==root){
        res_vect=new int[nsz*ar];
    }
    if(rank==root){
        for(int i=0;i<size;i++)szv[i]*=ar;
    }
    // cout << rank << " waits for gather" << endl;

    MPI_Gatherv(r.members,sz*ar,MPI_INT,res_vect,szv,displ,MPI_INT,root,MPI_COMM_WORLD);
    delete [] szv;
    delete [] arv;
    delete [] displ;

    // cout << rank << " through gather" << endl;

    if(rank!=root) {
        return *new relation();
    }
    else{
        return *new relation(res_vect,ar,nsz);
    }
}
relation& relation::distribute_itf(relation& r,hash_itf ith){
    int rank,size;int root=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    //get ar
    int ar=r.arity;
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
    if(ar==0)return *new relation();
    //get send lens and send displ
    int* slen=new int[size];
    for(int i=0;i<size;i++)slen[i]=0;
    for(int i=0;i<r.memsize;i++)slen[ith.get_value(r.get(i))]+=ar;
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
    int* data=new int[r.memsize*ar];
    for(int i=0;i<r.memsize;i++){
        int curid=ith.get_value(r.get(i));
        int curpos=ipos[curid]+sdis[curid];
        ipos[curid]+=ar;
        for(int j=0;j<ar;j++){
            data[curpos+j]=r.get(i)[j];
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

    return *new relation(collect,ar,ttl/ar);

}
relation& relation::distribute_loc(relation& r,hash_itf ith){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int ar=r.arity;

    vector<int> idset;

    for(int i=0;i<r.memsize;i++){
        int v=ith.get_value(r.get(i));
        if(v==rank)idset.push_back(i);
    }

    int nsz=idset.size();

    int* newmems=new int[nsz*ar];
    for(int i=0;i<nsz;i++){
        for(int j=0;j<ar;j++){
            newmems[i*ar+j]=r.get(idset[i])[j];
        }
    }

    return *new relation(newmems,ar,nsz);
}
relation& relation::distribute_hc(relation& r,int root,hash_hc cbh){

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);


    int l=r.memsize;
    MPI_Bcast(&l,1,MPI_INT,root,MPI_COMM_WORLD);
    int ar=r.arity;
    MPI_Bcast(&ar,1,MPI_INT,root,MPI_COMM_WORLD);

    //single source

    int* lenv=new int[size];
    for(int i=0;i<size;i++)lenv[i]=0;
    int* thresv=new int[size];

    if(rank==root){
        for(int i=0;i<r.memsize;i++){
            vector<int> cur=cbh.get_value(r.get(i));
            // print_array1d(cur);
            for(int j=0;j<cur.size();j++){
                lenv[cur[j]]+=ar;
            }
        }
        thresv[0]=0;
        for(int i=1;i<size;i++)thresv[i]=thresv[i-1]+lenv[i-1];
    }

    int len;
    MPI_Scatter(lenv,1,MPI_INT,&len,1,MPI_INT,root,MPI_COMM_WORLD);

    int* pos=new int[size];
    for(int i=0;i<size;i++)pos[i]=0;

    int* data=NULL;
    if(rank==root){

        int ttl=thresv[size-1]+lenv[size-1];
        data=new int[ttl];
        for(int i=0;i<r.memsize;i++){
            vector<int> cur=cbh.get_value(r.get(i));
            for(int j=0;j<cur.size();j++){
                int curid=cur[j];
                int curpos=thresv[curid]+pos[curid];
                pos[curid]+=ar;
                for(int k=0;k<ar;k++){
                    data[curpos+k]=r.get(i)[k];
                }
            }
        }

    }
    int* collect=new int[len];
    MPI_Scatterv(data,lenv,thresv,MPI_INT,collect,len,MPI_INT,root,MPI_COMM_WORLD);

    delete [] lenv;
    delete [] thresv;
    delete [] pos;
    delete [] data;

    return *new relation(collect,ar,len/ar);
}
relation& relation::distribute_hc_loc(relation& r,hash_hc cbh){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int ar=r.arity;

    vector<int> idset;

    for(int i=0;i<r.memsize;i++){
        vector<int> cur=cbh.get_value(r.get(i));
        for(int j=0;j<cur.size();j++){
            if(cur[j]==rank){
                idset.push_back(i);
                break;
            }
        }
    }

    int nsz=idset.size();

    int* newmems=new int[nsz*ar];
    for(int i=0;i<nsz;i++){
        for(int j=0;j<ar;j++){
            newmems[i*ar+j]=r.get(idset[i])[j];
        }
    }

    return *new relation(newmems,ar,nsz);
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

//core functions
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
    vector<vector<int> > idset;
    int p1=0,p2=0,l=comm.size();
    while(p1<f1.memsize&&p2<f2.memsize){
        int stat = cmpf(f1.get(p1),f2.get(p2),l,perms[0],perms[1]);
        if(stat > 0) p2++;
        if(stat < 0) p1++;
        if(stat == 0){
            int l1=0,l2=0;
            while(p1+l1<f1.memsize
                    &&cmpf(f1.get(p1),f1.get(p1+l1),l,perms[0],perms[0])==0)l1++;
            while(p2+l2<f2.memsize
                    &&cmpf(f2.get(p2),f2.get(p2+l2),l,perms[1],perms[1])==0)l2++;
            for(int i=0;i<l1;i++){
                for(int j=0;j<l2;j++){
                    vector<int> newmem(2);
                    newmem[0]=p1+i;
                    newmem[1]=p2+j;
                    idset.push_back(newmem);
                }
            }
            p1+=l1;
            p2+=l2;
        }
    }

    int newar=f1.arity+f2.arity-comm.size();
    int newsize=idset.size();
    int* newmems=new int[newar*newsize];
    for(int i=0;i<newsize;i++){
        int p1=idset[i][0];
        int p2=idset[i][1];
        for(int j=0;j<f1.arity;j++){
            newmems[i*newar+j]=f1.get(p1)[j];
        }
        int pos=f1.arity;
        for(int j=0;j<f2.arity;j++){
            if(is_commun[pat2.vars[j]])continue;
            newmems[i*newar+pos]=f2.get(p2)[j];
            pos++;
        }
    }
    return *new relation(newmems,newar,newsize);
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
relation& relation::join_mpi(relation& r1,relation& r2,
        pattern _pat1, pattern _pat2){
    int id, size; int root = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    pair<relation,relation::pattern> p1 = r1.filter(_pat1),p2 = r2.filter(_pat2);
    relation &f1 = p1.first, &f2 = p2.first;
    pattern pat1 = p1.second, pat2 = p2.second;
    vector<string> comm = pattern::find_comm(pat1,pat2);
    vector<int> id1,id2;
    int pos1=pat1.positions[comm[0]][0];
    int pos2=pat2.positions[comm[0]][0];

    for (int i=0; i<f1.memsize;i++){
        if (f1.get(i)[pos1]%size==id){
            id1.push_back(i);
        }
    }
    for (int i=0; i<f2.memsize;i++){
        if (f2.get(i)[pos2]%size==id){
            id2.push_back(i);
        }
    }
    int l1=id1.size(),l2=id2.size();
    int* mem1=new int[l1*f1.arity];
    for(int i=0;i<l1;i++){
        for(int j=0;j<f1.arity;j++){
            mem1[i*f1.arity+j]=f1.get(id1[i])[j];
        }
    }
    int* mem2=new int[l2*f2.arity];
    for(int i=0;i<l2;i++){
        for(int j=0;j<f2.arity;j++){
            mem2[i*f2.arity+j]=f2.get(id2[i])[j];
        }
    }
    relation &tmp1=*new relation(mem1,f1.arity,l1),&tmp2=*new relation(mem2,f2.arity,l2);
    relation &local= join(tmp1,tmp2,pat1,pat2);

    int len=local.memsize*local.arity;
    int* lenv=new int[size];
    MPI_Gather(&len,1,MPI_INT,lenv,1,MPI_INT,root,MPI_COMM_WORLD);
    int* disv=new int[size];
    disv[0]=0;for(int i=1;i<size;i++)disv[i]=disv[i-1]+lenv[i-1];
    int resl=disv[size-1]+lenv[size-1];
    int* newmems=new int[resl];
    MPI_Gatherv(local.members,len,MPI_INT,newmems,lenv,disv,MPI_INT,root,MPI_COMM_WORLD);
    delete [] lenv;
    delete [] disv;
    if(id!=root)return *new relation();
    else{
        int ar=f1.arity+f2.arity-comm.size();
        return *new relation(newmems,ar,resl/ar);
    }
}
relation& relation::join_itf(vector<relation>& rs,vector<pattern>& pats,int root){

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    bool quick_return=false;
    if(rank==root&&rs.size()==0){
        quick_return=true;
        MPI_Bcast(&quick_return,1,MPI_INT,root,MPI_COMM_WORLD);
    }
    if(quick_return)return *new relation();

    double ft=0,dt=0,jt=0,mt=0;

    #ifdef FLTR

    const clock_t s0=clock();
    pair<relation,relation::pattern> p=rs[0].filter(pats[0]);
    const clock_t t0=clock();
    ft+=(t0-s0)/(double)CLOCKS_PER_SEC;
    relation* currel=&p.first;
    pattern curpat=p.second;

    #else

    relation* currel=&rs[0];
    pattern curpat=pats[0];

    #endif

    hash_itf ith1(size),ith2(size);

    for(int i=1;i<rs.size();i++){

        #ifdef FLTR

        const clock_t s0=clock();
        pair<relation,relation::pattern> p=rs[i].filter(pats[i]);
        const clock_t t0=clock();
        ft+=(t0-s0)/(double)CLOCKS_PER_SEC;
        relation& nextrel=p.first;
        pattern nextpat=p.second;

        #else

        relation& nextrel=rs[i];
        pattern nextpat=pats[i];

        #endif

        vector<string> comm=pattern::find_comm(curpat,nextpat);
        int pos1=curpat.positions[comm[0]][0];
        int pos2=nextpat.positions[comm[0]][0];

        ith1.set_pos(pos1);
        ith2.set_pos(pos2);

        const clock_t s1=clock();
        relation& currel_local=distribute_itf(*currel,ith1);
        relation& next_local=distribute_loc(nextrel,ith2);
        const clock_t t1=clock();
        dt+=(t1-s1)/(double)CLOCKS_PER_SEC;

        const clock_t s2=clock();
        currel=&join(currel_local,next_local,curpat,nextpat);
        const clock_t t2=clock();
        jt+=(t2-s2)/(double)CLOCKS_PER_SEC;

        curpat=relation::pattern::join(curpat,nextpat);

        ith1.evolve();
        ith2.evolve();
    }
    const clock_t s3=clock();
    relation& res=merge(*currel,root);
    const clock_t t3=clock();
    mt+=(t3-s3)/(double)CLOCKS_PER_SEC;

    #ifdef FLTR

    if(rank==root)cout << "\tfltr:\t" << ft << endl;

    #endif

    if(rank==root)cout << "\tdtrb:\t" << dt << endl;
    if(rank==root)cout << "\tjoin:\t" << jt << endl;
    if(rank==root)cout << "\tmrge:\t" << mt << endl;

    return res;
}
relation& relation::join_hc(vector<relation>&_rs,vector<pattern>&_pats,int root){

    if(_rs.size()==0)return *new relation();

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    #ifdef FLTR

    vector<relation> rs;
    vector<pattern> pats;
    const clock_t s0=clock();
    for(int i=0;i<_rs.size();i++){
        pair<relation,pattern> p=_rs[i].filter(_pats[i]);
        rs.push_back(p.first);
        pats.push_back(p.second);
    }
    const clock_t t0=clock();
    if(rank==root)
        cout << "\tfltr:\t" << (t0-s0)/(double)CLOCKS_PER_SEC << endl;

    #else

    vector<relation>& rs=_rs;
    vector<pattern>& pats=_pats;

    #endif

    vector<string> comms=pattern::find_comm(pats);
    vector<int> ms=get_distribution(size,comms.size());

    vector<hash_hc> cbh;
    for(int i=0;i<pats.size();i++)
        cbh.push_back(*new hash_hc(ms,comms,pats[i]));

    vector<relation>& rs_loc=*new vector<relation>();

    const clock_t s1=clock();
    // for(int i=0;i<rs.size();i++)rs_loc.push_back(distribute_hc(rs[i],root,cbh[i]));
    for(int i=0;i<rs.size();i++)rs_loc.push_back(distribute_hc_loc(rs[i],cbh[i]));
    const clock_t t1=clock();

    if(rank==root)
        cout << "\tdtrb:\t" << (t1-s1)/(double)CLOCKS_PER_SEC << endl;

    const clock_t s2=clock();
    relation& res_loc=join(rs_loc,pats);
    const clock_t t2=clock();
    if(rank==root)
        cout << "\tjoin:\t" << (t2-s2)/(double)CLOCKS_PER_SEC << endl;

    const clock_t s3=clock();
    relation& res=merge(res_loc,root);
    const clock_t t3=clock();
    if(rank==root)
        cout << "\tmrge:\t" << (t3-s3)/(double)CLOCKS_PER_SEC << endl;

    return res;
}

//public interface function
relation& relation::join(relation& r1,relation& r2,string patstr1,string patstr2){
    return relation::join(r1,r2,*new pattern(patstr1),*new pattern(patstr2));
}
relation& relation::join(vector<relation>& rs,vector<string>& patstrs){
    if(rs.size()==1)return *new relation();
    vector<pattern> pats;
    for(int i=0;i<patstrs.size();i++){
        pats.push_back(pattern(patstrs[i]));
    }
    return join(rs,pats);
}
relation& relation::join_mpi(relation& r1,relation& r2,string patstr1,string patstr2){
    return join_mpi(r1,r2,pattern(patstr1),pattern(patstr2));
}
relation& relation::join_itf(vector<relation>&rs,vector<string>& strs){
    vector<pattern> pats;
    for(int i=0;i<strs.size();i++)
        pats.push_back(pattern(strs[i]));
    return relation::join_itf(rs,pats,0);
}
relation& relation::join_hc(vector<relation>&rs,vector<string>& patstrs){
    vector<pattern> pats;
    for(int i=0;i<patstrs.size();i++){
        pats.push_back(pattern(patstrs[i]));
    }
    return join_hc(rs,pats,0);
}

//output function
void relation::save(string filename){
    cout << "saving to " << filename << "..." << endl;
    ofstream file;
    file.open(("../output/"+filename).c_str());
    if(file.is_open()){
        file << arity << " " << memsize << "\n";
        for(int i=0;i<memsize;i++){
            for(int j=0;j<arity;j++){
                file << (members[i*arity+j]);
                if(j+1<arity) file << " ";
                else file << "\n";
            }
        }
        file.close();
    }
    else cout << "Saving failure." << endl;
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
ostream& operator<<(ostream &out, relation &rel){
    out << "#==begin==#" << endl;
    out << "ar:" << rel.arity  << " msz:" << rel.memsize << endl;
    for(int i=0;i<rel.memsize;i++){
        for(int j=0;j<rel.arity;j++){
            out << rel.get(i)[j] << ' ';
        }
        out << endl;
    }
    out << "ar:" << rel.arity << " msz:" << rel.memsize <<endl;
    out << "#==end==#" << endl;
    return out;
}
