#ifndef MEMBERDUMPER_HH_
#define MEMBERDUMPER_HH_

#include <memory>
#include <string>
#include <iostream>

template<class T>
class AbsMemberDumper
{
public:
    inline AbsMemberDumper(const std::string& i_label)
        : label_(i_label) {}

    inline virtual ~AbsMemberDumper() {}

    inline const std::string& label() const {return label_;}

    virtual void dump(const T& obj, std::ostream& os) const = 0;

private:
    std::string label_;
};

template<class T, class MemFcn>
class ConcreteMemberDumper : public AbsMemberDumper<T>
{
public:
    inline ConcreteMemberDumper(const std::string& i_label, MemFcn fcn)
        : AbsMemberDumper<T>(i_label), fcn_(fcn) {}

    inline virtual ~ConcreteMemberDumper() override {}

    virtual void dump(const T& obj, std::ostream& os) const override
        {os << fcn_(obj);}

private:
    MemFcn fcn_;
};

template<class T, class MemFcn>
std::shared_ptr<AbsMemberDumper<T> > make_MemberDumper(
    const std::string& i_label, MemFcn fcn)
{
    return std::shared_ptr<AbsMemberDumper<T> >(
        new ConcreteMemberDumper<T, MemFcn>(i_label, fcn));
}

template<class Coll>
void printMemberDumperLabels(const Coll& coll, std::ostream& os)
{
    const unsigned nDumpers = coll.size();
    for (unsigned i=0; i<nDumpers; ++i)
    {
        if (i) {os << ' ';}
        os << coll[i]->label();
    }
}

template<class Coll, class T>
void printMemberDumperValues(const Coll& coll, const T& obj, std::ostream& os)
{
    const unsigned nDumpers = coll.size();
    for (unsigned i=0; i<nDumpers; ++i)
    {
        if (i) {os << ' ';}
        coll[i]->dump(obj, os);
    }
}

#endif // MEMBERDUMPER_HH_
