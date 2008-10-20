#include "utl/defs.hh"
#include "propertiedClass.hpp"
#include "propertyValue.hpp"
#include <iostream>
#include <string>
using namespace std;

class NamedObj
{
public:
    NamedObj( string theName )
            :
            name( theName )
    {}

    string getName() const
    {
        return name;
    }

private:
    string name;
};

class Object : public NamedObj, public PropertiedClass
{
public:
    Object( string name )
            :
            NamedObj( name ),
            PropertiedClass()
    {}
};

int main()
{

    Object object( "Object 1" );
    object.createProperty( "Weight" );
    object.createProperty( "Mass" );

    object["Weight"] = Value( 220.0f );
    object["Height"] = Value( 72.0f );

    cout << object.getName() << endl;
    Value val = object["Weight"];
    cout << val << endl;
    object["Weight"] = object["Weight"].valueAsReal() + 10.0f;
    cout << object["Weight"] << endl;
    return 0;

}
