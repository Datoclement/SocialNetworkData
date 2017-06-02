#include <iostream>
#include <vector>
#include <string>
#include <map>

typedef int (hashtype)(int,int);

using namespace std;

class relation {

private:

    //arity, numbers of subsets that satisfy this relation
    int arity,memsize;

    //list of such subsets
    // vector<vector<int> > members;
    int* members;

    //private comparator used in sorting function
    struct cmp{
        vector<int> permutation;
        cmp(int arity);
        cmp(vector<int> _perm);
        bool operator()(int* a,int* b)const;
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

        //constructor from component
        pattern(vector<string> newvars,map<string,vector<int> > newpositions);

        //output stream
        friend ostream& operator<<(ostream& out,const pattern& pat);

        /**
            find commun variables in two patterns (used in join operation)
            @param  pat1, pat2: patterns
            @return           : list of commun variables
        */
        static vector<string> find_comm(pattern pat1,pattern pat2);

        /**
            find all commun variables that exist in more than 1 patterns
        */
        static vector<string> find_comm(vector<pattern> pats);

        /**
            find the corresponding permutation that is appropriate for the
            the join operation between pat1 and pat2
            @return: [perm for pat1; perm for pat2]
        */
        static vector<vector<int> > find_perm(pattern pat1,pattern pat2);

        /**
            find the result of pattern join
        */
        static pattern join(pattern pat1,pattern pat2);
    };

    /**
        a class that represents the complex hash function used in hypercube join
    */
    class hash_hc{

    private:

        /**
            mod number for each variable
            ordered according to the order in comm vector
        */
        vector<int> ms;

        /**
            the commun variables positions in corresponding pattern
            if it does not exist in the pattern, worth -1
            ordered by comm vector
        */
        vector<int> varspos;

        /**
            the corresponding group of processors when the hash value is given
        */
        vector<vector<int> >values;

    public:

        /**
            do three things:
            1. note down ms
            2. fill the table varspos in looking into comm and pat
            3. compute all groups of processors for each hash value
            @param ms  : mod numbers
            @param comm: commun variables of the group of studied pattern
            @param pat : the pattern for which we are constructing this hash class
        */
        hash_hc(vector<int>ms,vector<string>&comm,relation::pattern pat);

        /**
            compute the hash value of v, then given a list of processors
            that should receive it.
        */
        vector<int>& get_value(int* v);
    };

    /**
        hash class for instant-transfer algorithm
    */
    class hash_itf{

    private:
        int sz,r,pos;

    public:
        hash_itf(int _sz);

        void set_pos(int _pos);

        void evolve();

        int get_value(int* key);
    };

    //contructor with known components
    relation(int* _mem,int arity,int memsize);

    /**
        find relations that comply to the given pattern
        @param  pat: a given pattern
        @return    : a new relation object that contains all the relations that apply
    */
    pair<relation,pattern> filter(pattern pat);

    /**
        compare two tuples of different length with given permutations
        @param  l: compare first l elements after permutation
        @return  : 1 if a is lexicographically greater, -1 if less, 0 if equal
    */
    static int cmpf(int* a,int* b,int l,vector<int>& permutation_a,vector<int>& permutation_b);

    /**
        join two relations according to the given patterns
        @param  r1,r2    : relations
        @param  pat1,pat2: patterns
        @return          : a new relation object that contains the result
    */
    static relation& join(relation& r1, relation& r2, pattern pat1,pattern pat2);

    static relation& join(vector<relation>& rs,vector<string>&strs,
            relation&(*join_func)(vector<relation>&,vector<pattern>&,int));

    /**
        join several relations according to the given patterns
        implemented with a loop re-using simple 2-join
    */
    static relation& join_seq(vector<relation>& rs,vector<pattern>& pats,int root);

    /**
        join relations according to the given patterns
        implemented with naive algorithm with mpi
    */
    static relation& join_mpi(vector<relation>& rs,vector<pattern>& pats,int root);

    /**
        join several relations
        implemented with instant-transfer join algorithm
    */
    static relation& join_itf(vector<relation>& rs, vector<pattern>& pats, int root);

    /**
        join several relations
        implemented with hypercube join algorithm
    */
    static relation& join_hc(vector<relation>& rs, vector<pattern>& pats, int root);

    /**
        merge result from different machines
        @param  r   : local result
        @param  dest: destination of merge
        @return     : empty for non destination, merged result for destination
    */
    static relation& merge(relation& r,int dest);

    static relation& distribute_mpi(relation& r,int root,hash_itf ith);
    /**
        distribute intermedium result to corresponding processors
        essence of instant-transfer algorithm
        using MPI_Alltoall
        @param pos: pivot position that is taken into the hash function
    */
    static relation& distribute_itf(relation& r,hash_itf ith);

    /**
        distribute the raw data to corresponding processors
        essence of hypercube algorithm
        using MPI_Scatterv
        @param root: the root machine that has access to the raw data
        @param cbh : hash class (tuple->int)
    */
    static relation& distribute_hc(relation& r,hash_hc cbh);

    static relation& distribute_loc(relation& r,hash_itf ith);

public:

    //empty constructor
    relation();

    /**
        constructor from a given file
        @param str: the filename in the data directory
    */
    relation(string str);

    //output stream functions
    friend ostream& operator<<(ostream& out, relation& rel);
    friend ostream& operator<<(ostream& out,relation::pattern& pat);

    int* get(int id);

    //sort the members by built-in comparator of vector class
    void sort();

    //sort the members with a given permutation
    void sort(const vector<int>& permutation);

    void sort(cmp f);

    //seeding of the random number generator
    void random_seed(int seed);

    //return a random member
    // vector<int> random();
    int* random();

    //save the relation into the given path
    void save(string filename);

    void free();

    //join several relations with strings that correspond to a pattern
    static relation& join_seq(vector<relation>& rs,vector<string>& patstrs);

    //join several relations with naive algo
    static relation& join_mpi(vector<relation>& rs,vector<string>& patstrs);

    //join several relations with itf algo
    static relation& join_itf(vector<relation>& rs,vector<string>& patstrs);

    //join several relation with hc algo
    static relation& join_hc(vector<relation>& rs, vector<string>& patstrs);

};
