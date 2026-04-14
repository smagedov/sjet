#ifndef MINIMALEVENT_HH_
#define MINIMALEVENT_HH_

#include <limits>
#include <stdexcept>

// Event structures playing well with the analysis framework
// should have the methods described below:
//
// Default constructor.
//
// Static method "version".
//
// Methods for dealing with event numbers: "setNumber", "number",
// "isNumberValid".
//
// Method "clear" for invalidating the event contents. It is expected
// that the event structure will be created only once before entering
// the event loop and then the "clear" and "setNumber" methods will
// be called inside the loop before everything else happens.
//
struct MinimalEvent
{
    inline MinimalEvent() {clear();}

    // Data contained in this event (just an example)
    double randomDatum;
    bool randomDatumReady;

    // Handle event numbers
    inline void setNumber(const unsigned long n)
    {
        if (n == invalidEventNumber_)
            throw std::invalid_argument("In MinimalEvent::setNumber : "
                                        "invalid event number");
        number_ = n;
    }
    inline unsigned long number() const {return number_;}

    // Check if the event number has been set
    inline bool isNumberValid() const
        {return number_ != invalidEventNumber_;}

    // Invalidate event number and various data sections
    void clear()
    {
        number_ = invalidEventNumber_;
        randomDatumReady = false;
    }

    // Bump up the version number if you change
    // event contents
    static inline unsigned version() {return 1;}

private:
    static constexpr unsigned long invalidEventNumber_ =
        std::numeric_limits<unsigned long>::max();
    unsigned long number_;
};

#endif // MINIMALEVENT_HH_
