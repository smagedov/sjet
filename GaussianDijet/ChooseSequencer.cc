#include <cassert>
#include <limits>

#include "ChooseSequencer.hh"

ChooseSequencer::ChooseSequencer(const unsigned m_in, const unsigned n_in)
    : state_(n_in, 0U), m_(m_in), n_(n_in), valid_(true)
{
    assert(m_ <= n_);
    assert(n_ < std::numeric_limits<unsigned>::max());

    for (unsigned i=0; i<m_; ++i)
        state_[i] = i+1U;
}

void ChooseSequencer::getChoice(unsigned* choice, const unsigned sz) const
{
    assert(valid_);
    assert(sz >= m_);

    unsigned ind = 0U;
    for (unsigned i=0U; i<n_; ++i)
        if (state_[i])
            choice[ind++] = i;
}

ChooseSequencer& ChooseSequencer::operator++()
{
    if (valid_)
    {
        if (m_ == 0 || m_ == n_)
            valid_ = false;
        else
        {
            // Find the index of the first 0 from the right
            unsigned fInd = n_ - 1U;
            while (state_[fInd])
                --fInd;
            if (fInd < n_ - m_)
                valid_ = false;
            else
            {
                // Find the index of the first non-0 below that 0
                unsigned non0 = fInd - 1U;
                while (!state_[non0])
                    --non0;

                // Move that non-0 element to the right and
                // renumber the elements above the moved one
                unsigned elem = state_[non0];
                state_[non0++] = 0U;
                for (; elem<=m_; ++elem)
                    state_[non0++] = elem;
                for (; non0<n_; ++non0)
                    state_[non0] = 0U;
            }
        }
    }
    return *this;
}
