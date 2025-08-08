#ifndef SJET_ABSNODEVISITOR_HH_
#define SJET_ABSNODEVISITOR_HH_

namespace sjet {
    template<class Particle>
    struct AbsNodeVisitor
    {
        typedef Particle particle_type;

        inline virtual ~AbsNodeVisitor() {}

        virtual void visit(unsigned nodeIndex, const Particle& p) = 0;
    };
}

#endif // SJET_ABSNODEVISITOR_HH_
