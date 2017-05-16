#include <iostream>
#include <vector>
#include <string>
#include <map>

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


        pattern(vector<string> newvars,map<string,vector<int> > newpositions);
    public:
        int size;

        //variables
        vector<string> vars;

        //map from variable name to its postion
        map<string,vector<int> > positions;

        //constructor from a string
        pattern(string pat);

        //output stream
        friend ostream& operator<<(ostream& out,const pattern& pat);

        /**
            find commun variables in two patterns (used in join operation)
            @param  pat1, pat2: patterns
            @return           : list of commun variables
        */
        static vector<string> find_comm(pattern pat1,pattern pat2);

        /**
            find the corresponding permutation that is appropriate for the
            the join operation between pat1 and pat2
            @return: [perm for pat1; perm for pat2]
        */
        static vector<vector<int> > find_perm(pattern pat1,pattern pat2);

        
        static pattern join(pattern pat1,pattern pat2);
    };

    //empty constructor
    relation();

    //contructor with a given set of relations of the same arity
    relation(vector<vector<int> > _mem);

    /**
        find relations that comply to the given pattern
        @param  pat: a given pattern
        @return    : a new relation object that contains all the relations that apply
    */
    relation filter(pattern pat);

    /**
        join two relations according to the given patterns
        @param  r1,r2    : relations
        @param  pat1,pat2: patterns
        @return          : a new relation object that contains the result
    */
    static relation join(relation r1, relation r2, pattern pat1,pattern pat2);

public:

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
    static relation join(relation r1, relation r2, string patstr1, string patstr2);

    //join several relations with strings that correspond to a pattern
    static relation join(vector<relation> rs,vector<string> patstrs);
    // static relation multijoin(vector<relation> rs,vector<string> patstrs);
};
