#include <iostream>
#include <vector>
#include <string>
#include <map>

typedef int (hashtype)(int,int);

using namespace std;

class relation {

private:

    //identifier
    string id;

    //arity, numbers of subsets that satisfy this relation
    int arity,memsize;

    //list of such subsets
    vector<vector<int> > members;

    //private comparator used in sorting function
    struct cmp{
        vector<int> permutation;
        cmp(int arity);
        cmp(vector<int> _perm);
        bool operator()(vector<int> a,vector<int> b)const;
    };

    /**
        a class that represents the relationship between variables of relations
        and their positions
    */
    class pattern{
    private:

        //only a~z and non-beginning 1~9 is valid character in a variable name
        bool is_valid_char(char c);

    public:

        int size;

        //variables
        vector<string> vars;

        //map from variable name to its postion
        map<string,vector<int> > positions;

        //constructor from a string
        pattern(string pat);

        pattern(vector<string> newvars,map<string,vector<int> > newpositions);

        //output stream
        friend ostream& operator<<(ostream& out,const pattern& pat);

        /**
            find commun variables in two patterns (used in join operation)
            @param  pat1, pat2: patterns
            @return           : list of commun variables
        */
        static vector<string> find_comm(pattern pat1,pattern pat2);

        static vector<string> find_comm(vector<pattern> pats);

        /**
            find the corresponding permutation that is appropriate for the
            the join operation between pat1 and pat2
            @return: [perm for pat1; perm for pat2]
        */
        static vector<vector<int> > find_perm(pattern pat1,pattern pat2);


        static pattern join(pattern pat1,pattern pat2);
    };

    class cubehash{

    private:
        vector<int> ms;

        vector<int> varspos;

        vector<vector<int> >values;

    public:
        cubehash(vector<int>ms,vector<string>&comm,relation::pattern pat);

        vector<int>& get_value(vector<int>& v);
    };

    //contructor with a given set of relations of the same arity
    relation(vector<vector<int> >& _mem);

    /**
        find relations that comply to the given pattern
        @param  pat: a given pattern
        @return    : a new relation object that contains all the relations that apply
    */
    pair<relation,pattern> filter(pattern pat);

    /**
        join two relations according to the given patterns
        @param  r1,r2    : relations
        @param  pat1,pat2: patterns
        @return          : a new relation object that contains the result
    */
    static relation& join(relation& r1, relation& r2, pattern pat1,pattern pat2);

    static relation& join(vector<relation>& rs,vector<pattern>& pats);

    /**
        join two relations according to the given patterns
    */
    static relation& join_mpi_hash(relation& r1, relation& r2, pattern p1, pattern p2, hashtype* hash);

    static relation& merge_vector_to_relation(vector<vector<int> >&mems, int dest);

    // static vector<vector<int> >& distribute(vector<vector<int> >&mems,int source,int pos);
    static vector<vector<int> >& distribute(vector<vector<int> >&mems,int pos);


    static vector<vector<int> >& distribute_cube(vector<vector<int> >&mems,int source,cubehash chb);

    static relation& join_mpi(vector<relation>& rs, vector<pattern>& pats, int root);

    static relation& join_hypercube(vector<relation>& rs, vector<pattern>& pats, int root);

public:

    //empty constructor
    relation();

    /**
        constructor from a given file
        @param str: the filename in the data directory
    */
    relation(string str);

    //output stream functions
    friend ostream& operator<<(ostream& out, const relation& rel);
    friend ostream& operator<<(ostream& out,relation::pattern& pat);

    //sort the members by built-in comparator of vector class
    void sort();

    //sort the members with a given permutation
    void sort(const vector<int>& permutation);

    //seeding of the random number generator
    void random_seed(int seed);

    //return a random member
    vector<int> random();

    //save the relation into the given path
    void save(string filename);

    //join two relations with two strings that correspond to a pattern
    static relation& join(relation& r1, relation& r2, string patstr1, string patstr2);

    //join two relations using mpi
    static relation& join_mpi(relation& r1, relation& r2, string patstr1, string patstr2);

    //join several relations with strings that correspond to a pattern
    static relation& join(vector<relation>& rs,vector<string>& patstrs);

    //join several relations using mpi with strings that correspond to a pattern
    static relation& join_mpi(vector<relation>& rs,vector<string>& patstrs);

    static relation& join_hypercube(vector<relation>& rs, vector<string>& pats);

};
