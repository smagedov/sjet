#ifndef SJET_NODEVISITORS_HH_
#define SJET_NODEVISITORS_HH_

#include <vector>
#include <algorithm>

#include "sjet/AbsNodeVisitor.hh"

namespace sjet {
    // Visitor that adds particle 4-vectors
    template<class Particle>
    class ParticleAdder : public AbsNodeVisitor<Particle>
    {
    public:
        inline ParticleAdder() {}
        inline virtual ~ParticleAdder() override {}

        inline virtual void visit(unsigned /* node */, const Particle& p) override
            {sum_ += p;}

        inline const Particle& result() const {return sum_;}

    private:
        Particle sum_;
    };

    // Visitor that adds particle scalar pt values
    template<class Particle>
    class ScalarPtAdder : public AbsNodeVisitor<Particle>
    {
    public:
        inline ScalarPtAdder() : sum_(0.0L) {}
        inline virtual ~ScalarPtAdder() override {}

        inline virtual void visit(unsigned /* node */, const Particle& p) override
            {sum_ += p.pt();}

        inline double result() const {return sum_;}

    private:
        long double sum_;
    };

    // Visitor that collects the indices of the visited nodes
    template<class Particle>
    class NodeCollector : public AbsNodeVisitor<Particle>
    {
    public:
        inline NodeCollector() {}
        inline virtual ~NodeCollector() override {}

        inline virtual void visit(const unsigned node, const Particle& /* p */) override
            {nodes_.push_back(node);}

        inline void sort()
            {std::sort(nodes_.begin(), nodes_.end());}

        inline const std::vector<unsigned>& result() const {return nodes_;}

    private:
        std::vector<unsigned> nodes_;
    };

    // Visitor that simply counts the visited nodes
    template<class Particle>
    class NodeCounter : public AbsNodeVisitor<Particle>
    {
    public:
        inline NodeCounter() : count_(0) {}
        inline virtual ~NodeCounter() override {}

        inline virtual void visit(unsigned /* node */, const Particle& /* p */) override
            {++count_;}

        inline unsigned result() const {return count_;}

    private:
        unsigned count_;
    };
}

#endif // SJET_NODEVISITORS_HH_
