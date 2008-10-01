class ExpressionEvaluator
{
protected:
    typedef ExpressionCompiler::Code Code;
    
    typedef void* Pointer;

    class VirtualMachine
    {
        union StackElement
        {
            Real    real;
            Pointer pointer;
            Integer integer;
        };
    public:
        VirtualMachine()
        {}

        ~VirtualMachine()
        {}

        const Real execute( const Code& code );
    };

    VirtualMachine virtualMachine;

    bool recompileFlag;
    
    


};
