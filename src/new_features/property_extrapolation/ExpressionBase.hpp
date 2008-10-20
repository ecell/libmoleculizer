#include "ExpressionCompiler.hpp"


class ExpressionBase
{
protected:
    typedef ExpressionCompiler::Code Code;

    typedef void* Pointer;

    class VirtualMachine
    {
        union StackElement
        {
            Real real;
            Pointer pointer;
            Integer integer;
        };

    public:
        VirtualMachine()
        {}

        ~VirtualMachine()
        {}

        Real execute( const Code& code );

        void setExpression( const std::string& value );
        const std::string& getExpression() const;



    protected:
        std::string expression;
        bool recompileFlag;

    private:

    };

};
