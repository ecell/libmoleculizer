#include

typedef flaot Real;

// This class is generic, like the Expression Evaluator.
// This class is meant to represent

class Extrapolatable : public ExpressionEvaluator
{
public:
    Extrapolatable();
    ~Extrapolatable();

    virtual Real getValue()
    {
        ExpressionCompiler compiler( this, propertyMap );

        Code compiledCode;
        compiledCode = compiler.compileExpression( this->expression );

        Real value  = virtualMachine.execute( compiledCode );

        return value;
    }

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
    ExpressionEvaluator

    bool recompileFlag;

};

const Real ExpressionProcessBase::VirtualMachine::execute( const Code& code )
{

#define FETCH_OPCODE()\
    reinterpret_cast<const ExpressionCompiler::InstructionHead* const>( pC )\
      ->getOpcode()

#define DECODE_INSTRUCTION( OPCODE )\
    typedef ExpressionCompiler::\
      Opcode2Instruction<ExpressionCompiler::OPCODE>::type CurrentInstruction;\
    const CurrentInstruction* const instruction\
      ( reinterpret_cast<const CurrentInstruction* const>( pC ) )


#define INCREMENT_PC( OPCODE )\
    pC += sizeof( ExpressionCompiler::\
                   Opcode2Instruction<ExpressionCompiler::OPCODE>::type );\

  //    std::cout << #OPCODE << std::endl;

    StackElement stack[100];
    //  stack[0].real = 0.0;
    StackElement* stackPtr( stack - 1 );

    const unsigned char* pC( &code[0] );

    while ( 1 )
    {

        Real bypass;

        switch ( FETCH_OPCODE() )
        {

#define SIMPLE_ARITHMETIC( OPCODE, OP )\
     ( stackPtr - 1)->real OP##= stackPtr->real;\
     INCREMENT_PC( OPCODE );\
     --stackPtr

            /*
                     const Real topValue( stackPtr->real );\
              INCREMENT_PC( OPCODE );\
              ( stackPtr - 1 )->real OP##= topValue;\
              --stackPtr;\
            */

        case ExpressionCompiler::ADD:
        {
            SIMPLE_ARITHMETIC( ADD, + );

            continue;
        }

        case ExpressionCompiler::SUB:
        {
            SIMPLE_ARITHMETIC( SUB, - );

            continue;
        }

        case ExpressionCompiler::MUL:
        {
            SIMPLE_ARITHMETIC( MUL, * );

            continue;
        }

        case ExpressionCompiler::DIV:
        {
            SIMPLE_ARITHMETIC( DIV, / );

            continue;
        }

#undef SIMPLE_ARITHMETIC

        case ExpressionCompiler::CALL_FUNC2:
        {
            DECODE_INSTRUCTION( CALL_FUNC2 );

            ( stackPtr - 1 )->real
            = ( instruction->getOperand() )(( stackPtr - 1 )->real,
                                            stackPtr->real );
            --stackPtr;

            INCREMENT_PC( CALL_FUNC2 );
            continue;
        }


        case ExpressionCompiler::CALL_FUNC1:
        {
            DECODE_INSTRUCTION( CALL_FUNC1 );

            stackPtr->real
            = ( instruction->getOperand() )( stackPtr->real );

            INCREMENT_PC( CALL_FUNC1 );
            continue;
        }

        case ExpressionCompiler::NEG:
        {
            stackPtr->real = - stackPtr->real;

            INCREMENT_PC( NEG );
            continue;
        }

#if 0
        case ExpressionCompiler::PUSH_INTEGER:
        {
            DECODE_INSTRUCTION( PUSH_INTEGER );

            ++stackPtr;
            stackPtr->integer = instruction->getOperand();

            INCREMENT_PC( PUSH_INTEGER );
            continue;
        }

        case ExpressionCompiler::PUSH_POINTER:
        {
            DECODE_INSTRUCTION( PUSH_POINTER );

            ++stackPtr;
            stackPtr->pointer = instruction->getOperand();

            INCREMENT_PC( PUSH_POINTER );
            continue;
        }

#endif // 0

        case ExpressionCompiler::PUSH_REAL:
        {
            DECODE_INSTRUCTION( PUSH_REAL );

            bypass = instruction->getOperand();

            INCREMENT_PC( PUSH_REAL );
            goto bypass_real;
        }

        case ExpressionCompiler::LOAD_REAL:
        {
            DECODE_INSTRUCTION( LOAD_REAL );

            bypass = *( instruction->getOperand() );

            INCREMENT_PC( LOAD_REAL );
            goto bypass_real;
        }

        case ExpressionCompiler::OBJECT_METHOD_REAL:
        {
            DECODE_INSTRUCTION( OBJECT_METHOD_REAL );

            bypass = ( instruction->getOperand() )();

            INCREMENT_PC( OBJECT_METHOD_REAL );
            goto bypass_real;
        }

        case ExpressionCompiler::OBJECT_METHOD_INTEGER:
        {
            DECODE_INSTRUCTION( OBJECT_METHOD_INTEGER );

            bypass = static_cast<Real>(( instruction->getOperand() )() );

            INCREMENT_PC( OBJECT_METHOD_INTEGER );
            goto bypass_real;
        }

        case ExpressionCompiler::RET:
        {
            return stackPtr->real;
        }

        default:
        {
            THROW_EXCEPTION( UnexpectedError, "Invalid instruction." );
        }

        }

#if defined( ENABLE_STACKOPS_FOLDING )

bypass_real:

        // Fetch next opcode, and if it is the target of of the stackops folding,
        // do it here.   If not (default case), start the next loop iteration.
        switch ( FETCH_OPCODE() )
        {
        case ExpressionCompiler::ADD:
        {
            stackPtr->real += bypass;

            INCREMENT_PC( ADD );
            break;
        }

        case ExpressionCompiler::SUB:
        {
            stackPtr->real -= bypass;

            INCREMENT_PC( SUB );
            break;
        }

        case ExpressionCompiler::MUL:
        {
            stackPtr->real *= bypass;

            INCREMENT_PC( MUL );
            break;
        }

        case ExpressionCompiler::DIV:
        {
            stackPtr->real /= bypass;

            INCREMENT_PC( DIV );
            break;
        }

        case ExpressionCompiler::CALL_FUNC2:
        {
            DECODE_INSTRUCTION( CALL_FUNC2 );

            stackPtr->real
            = ( instruction->getOperand() )( stackPtr->real, bypass );

            INCREMENT_PC( CALL_FUNC2 );
            break;
        }

        default:
        {
            // no need to do invalid instruction check here because
            // it will be done in the next cycle.

            ++stackPtr;
            stackPtr->real = bypass;

            break;
        }
        }

        continue;

#else /* defined( ENABLE_STACKOPS_FOLDING ) */

bypass_real:

        ++stackPtr;
        stackPtr->real = bypass;

        continue;

#endif /* defined( ENABLE_STACKOPS_FOLDING ) */

    }

#undef DECODE_INSTRUCTION
#undef FETCH_INSTRUCTION
#undef INCREMENT_PC

}

