#ifndef CHOOSESEQUENCER_HH_
#define CHOOSESEQUENCER_HH_

//===========================================================================
// ChooseSequencer.hh
//
// The class below will generate all combinations of m elements
// drawn without replacement from a set of n elements. Naturally,
// there are n!/(m! (n-m)!) such choices (that is, n choose m).
//
// The intended use is as follows:
//
// std::vector<unsigned> chosenElements(m);
//
// for (ChooseSequencer seq(m, n); seq.isValid(); ++seq)
// {
//     seq.getChoice(m ? &chosenElements[0] : nullptr, chosenElements.size());
//
//     ... At this point, the "chosenElements" vector will be filled ......
//     ... with m distinct numbers, each between 0 and n-1 (inclusive) ....
//     ... in such a way that the set of numbers will be different for ....
//     ... every iteration. These m numbers will be arranged in the .......
//     ... increasing order. Now, do whatever needs to be done with .......
//     ... the set chosen in this iteration. ..............................
// }
//
// I. Volobouev
// February 2026
//===========================================================================

#include <vector>

class ChooseSequencer
{
public:
    // Choosing m out of n. Must have m <= n.
    ChooseSequencer(unsigned m, unsigned n);

    inline unsigned m() const {return m_;}
    inline unsigned n() const {return n_;}
    inline bool isValid() const {return valid_;}

    // "choice" must point to an array with at least m elements
    void getChoice(unsigned* choice, unsigned choiceArraySize) const;

    ChooseSequencer& operator++();

private:
    std::vector<unsigned> state_;
    unsigned m_;
    unsigned n_;
    bool valid_;
};

#endif // CHOOSESEQUENCER_HH_
