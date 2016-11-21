#ifndef AMSimulation_CombinationFactory_h_
#define AMSimulation_CombinationFactory_h_

#include <vector>

namespace slhcl1tt {

class CombinationFactory {
  public:
    // Constructor
    CombinationFactory() : verbose_(true) {}

    // Destructor
    ~CombinationFactory() {}

    // Enum
    enum Flag { BAD=999999999 };

    // Arrange combinations
    // groups[i][j] is the j-th element in the i-th group
    // combinations[i][j] is the j-th element in the i-th combination
    template<typename T>
    std::vector<std::vector<T> > really_combine(const std::vector<std::vector<T> >& groups) {
        std::vector<T> combination;
        std::vector<std::vector<T> > combinations;

        const int ngroups = groups.size();
        std::vector<unsigned> indices(ngroups, 0);  // init to zeroes

        int i=0, j=0;
        while (true) {
            combination.clear();
            for (i=0; i<ngroups; ++i) {
                if (groups.at(i).size())
                    combination.push_back(groups.at(i).at(indices.at(i)));
                else  // empty group
                    combination.push_back(CombinationFactory::BAD);
            }
            combinations.push_back(combination);

            for (i=ngroups-1; i>=0; --i)
                if (groups.at(i).size())
                    if (indices.at(i) != groups.at(i).size() - 1)
                        break;  // take the last index that has not reached the end
            if (i == -1)  break;

            indices[i] += 1;  // increment that index
            for (j=i+1; j<ngroups; ++j)
                indices[j] = 0;  // set indices behind that index to zeroes
        }

        return combinations;
    }

    template<typename T>
    std::vector<std::vector<T> > combine(const std::vector<std::vector<T> >& groups, bool FiveOfSix) {
        if (!FiveOfSix)  return really_combine(groups);  // do not have to do 5/6 combinations in 6/6

        assert(groups.size() == 6);
        int cnt = 0;
        for (unsigned i=0; i<groups.size(); ++i) {
            cnt += (groups.at(i).size() > 0);
        }
        if (cnt != 6)  return really_combine(groups);  // not 6/6

        std::vector<std::vector<T> > combinations = really_combine(groups);  // 6/6
        std::vector<std::vector<T> > tmp_groups = groups;  // clone

        for (unsigned i=0; i<groups.size(); ++i) {
            tmp_groups.at(i).clear();
            const std::vector<std::vector<T> >& tmp_combinations = really_combine(tmp_groups);  // 5/6
            combinations.insert(combinations.end(), tmp_combinations.begin(), tmp_combinations.end());
            tmp_groups.at(i) = groups.at(i);
        }
        return combinations;
    }

    // Debug
    void print();


  private:
    // Member data
    int verbose_;
};

}

#endif
