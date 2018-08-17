
#include <vector>
#include <map>

extern "C"
{
#include "family.h"
}

namespace wgd
{
    //conditioning = "oneOrMore" if all families with one or more gene copies are included in the data,
    //conditioning = "twoOrMore" to condition on families having two of more genes,
    //conditioning = "oneInBothClades" if the data set was filtered to include only families with at least one gene copy
    //in each of the two main clades stemming from the root.
    //conditioning = "none" uses unconditional likelihoods.
    enum conditioning { none, oneOrMore, twoOrMore, oneInBothClades };

    struct wgd
    {
        std::vector<double> retain1;
        std::vector<double> retain2;
        std::vector<double> retain3;
    };

    struct Rate
    {
        double retention;
        double loss;
    };

    struct Doom
    {
        double doom;
        double doom_left;
        double doom_right;

        void set(bool left, double val)
        {
            if (left)
                doom_left = val;
            else
                doom_right = val;
        }

        void add(bool left, double val)
        {
            if (left)
                doom_left += val;
            else
                doom_right += val;
        }

        void update(bool left)
        {
            if (left)
                doom = doom_left;
            else
                doom = doom_right;
        }
    };
    struct phylo_data
    {
        //        enum Types { Unknown, BirthDeath, WholeGenomeDuplicate, WholeGenomeTriplicate, RootPrior };
        //        pCafeNode *Parent;
        //        pCafeNode *Child;
        double Time;
        //        double RetenRate;
        //        double LossRate;
        //        Types Type;
        std::string Species;
    };

    struct wgd_input
    {
#if 0
        std::map<pCafeNode, phylo_data> phyloMat;
        double para;
        double lower;
        double upper;
        wgd wgdTab;
#endif
    };


    double getLikGeneCount(std::vector<double> para, const wgd_input& b, pCafeFamily fam, int mMax, double geomProb, conditioning c);
    wgd_input processInput(pCafeTree pTree, double startingBDrate);
    int binomialCoeff(int n, int k);
    double calcProbSamePara(double logLam, double ti, int mStart, int nEnd);
    double calcProb(std::vector<double> logLamlogMu, double ti, int mStart, int nEnd);
    std::vector<std::vector<double>> getPt(std::vector<double> logLamlogMu, double ti, int nPos, int nFamily, bool isChildLeaf, std::vector<int> nLeafCount, bool isWgdNode, double wgdLR);
    void calcProbOneBranch(std::vector<pCafeNode> branchNode,
        std::vector<double> logLamlogMu,
        int nPos,
        int nFamily,
        int nLeaf,
        std::map<std::string, std::vector<int>> geneCountData,
        std::map<pCafeNode, Rate> wgdTab,
        std::map<pCafeNode, phylo_data> phyloMat,
        std::map<pCafeNode, Doom> MatDoomed,
        bool doomLeft,
        std::vector<std::vector<std::vector<double>>> Mat);
}