#include <boost/function.hpp>


template <typename T>
class callback
{
private:
    boost::function<void (T)> _f;

public:
    command()
    {}
    
    command(boost::function<void ()> f) 
        : 
        _f(f) 
    {}

    template <typename T> 
    void setFunction (T t) 
    {
        _f = t;
    }

    void execute()
    {
        if(!_f.empty())
            _f();
    }
};
