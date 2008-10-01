//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//       This file is part of the E-Cell System
//
//       Copyright (C) 1996-2008 Keio University
//       Copyright (C) 2005-2008 The Molecular Sciences Institute
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// E-Cell is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
//END_HEADER
//
// Original Authors:
//     Koichi Takahashi
//     Tatsuya Ishida
//
// Modifying Authors:
//     Nathan Addy, Molecular Sciences Institute
//


#ifndef __EXPRESSIONCOMPILER_HPP
#define __EXPRESSIONCOMPILER_HPP

#include <new>

#include <boost/spirit/core.hpp>
#include <boost/spirit/tree/ast.hpp>

#if SPIRIT_VERSION >= 0x1800
#define PARSER_CONTEXT parser_context<>
#else
#define PARSER_CONTEXT parser_context
#endif

#include "utl/AssocVector.h"

// #include "libecs/libecs.hpp"
// #include "libecs/Process.hpp"
// #include "libecs/MethodProxy.hpp"

using namespace boost::spirit;

// Pointers to function that match a string to a real func(real) pointer.
typedef
::Loki::AssocVector < String, Real( * )( Real ), std::less<const String> >
FunctionMap1;

typedef
::Loki::AssocVector < String, Real( * )( Real, Real ), std::less<const String> >
FunctionMap2;

typedef
::Loki::AssocVector<String, Real, std::less<const String> > ConstantMap;

typedef
::Loki::AssocVector< String, Real, std::less<const String> > PropertyMap;

class ExpressionCompiler
{
public:

    typedef std::vector<unsigned char> Code;

    typedef std::vector<char> CharVector;

    // possible operand types:
    typedef utl::Real Real;
    typedef Real* RealPtr;
    typedef utl::Integer Integer;
    typedef void* Pointer;

//     typedef Entity*  ( libecs::VariableReference::* VariableReferenceEntityMethodPtr )() const;
//     typedef Entity*  ( libecs::Process::* ProcessMethodPtr )() const;
//    typedef const Real ( libecs::Entity::* EntityMethodPtr )() const;

    typedef Real ( *RealFunc0 )();
    typedef Real ( *RealFunc1 )( Real );
    typedef Real ( *RealFunc2 )( Real, Real );


    typedef ObjectMethodProxy<Real> RealObjectMethodProxy;
    typedef ObjectMethodProxy<Integer> IntegerObjectMethodProxy;
    //  VariableReferenceMethodProxy;

    typedef void ( *InstructionAppender )( Code& );


    enum Opcode// the order of items is optimized. don't change.
    {
        ADD = 0  // no arg
        , SUB    // no arg
        , MUL    // no arg
        , DIV    // no arg
        , CALL_FUNC2 // RealFunc2
        // Those instructions above are candidates of stack operations folding
        // in the stack machine, and grouped here to optimize the switch().

        , CALL_FUNC1 // RealFunc1
        //, CALL_FUNC0 // RealFunc0
        , NEG    // no arg
        //      , PUSH_POINTER // Pointer
        , PUSH_REAL   // Real
        , LOAD_REAL  // Real*
        , OBJECT_METHOD_REAL //VariableReferencePtr, VariableReferenceMethodPtr
        , OBJECT_METHOD_INTEGER // VariableReferenceIntegerMethodPtr
        , RET   // no arg
        , END = RET
    };

    class NoOperand
        {}; // Null type.

    template <Opcode OPCODE>
    class Opcode2Operand
    {
    public:
        typedef NoOperand type;
    };


    template <Opcode OPCODE>
    class Opcode2Instruction;

    class InstructionHead
    {
    public:

        InstructionHead( Opcode opcode )
                :
                opcode( opcode )
        {
            ; // do nothing
        }

        const Opcode getOpcode() const
        {
            return this->opcode;
        }

    private:

        const Opcode  opcode;

    };

    template < class OPERAND >
    class InstructionBase
                :
                public InstructionHead
    {
    public:

        InstructionBase( Opcode opcode, const OPERAND& operand )
                :
                InstructionHead( opcode ),
                operand( operand )
        {
            ; // do nothing
        }

        const OPERAND& getOperand() const
        {
            return this->operand;
        }

    private:

        InstructionBase( Opcode );

    protected:

        const OPERAND operand;
    };


    /**
       Instruction Class
    */

    template < Opcode OPCODE >
    class Instruction
                :
                public InstructionBase<typename Opcode2Operand<OPCODE>::type >
    {

    public:

        typedef typename Opcode2Operand<OPCODE>::type Operand;

        Instruction( const Operand& operand )
                :
                InstructionBase<Operand>( OPCODE, operand )
        {
            ; // do nothing
        }

        Instruction()
                :
                InstructionBase<Operand>( OPCODE )
        {
            ; // do nothing
        }

    };


    template <class CLASS, typename RESULT>
    struct ObjectMethodOperand
    {
        //typedef boost::mem_fn< RESULT, CLASS > MethodType;
        typedef RESULT ( CLASS::* MethodPtr )( void ) const;

        const CLASS* operand1;
        MethodPtr operand2;
    };

    typedef ObjectMethodOperand<Process, Real> ProcessMethod;
    typedef ObjectMethodOperand<Entity, Real>  EntityMethod;


private:

    class CompileGrammar
                :
                public grammar<CompileGrammar>
    {
    public:
        enum GrammarType
        {
            GROUP = 1,
            INTEGER,
            FLOAT,
            NEGATIVE,
            EXPONENT,
            FACTOR,
            POWER,
            TERM,
            EXPRESSION,
            VARIABLE,
            CALL_FUNC,
            ENTITY_FUNC,
            ENTITY_PROPERTY,
            IDENTIFIER,
            CONSTANT,
        };

        template <typename ScannerT>
        struct definition
        {
#define leafNode( str ) leaf_node_d[lexeme_d[str]]
#define rootNode( str ) root_node_d[lexeme_d[str]]

            definition( CompileGrammar const& /*self*/ )
            {
                integer     =   leafNode( + digit_p );
                floating    =   leafNode( + digit_p >> ch_p( '.' ) >> + digit_p );

                exponent    =   ( floating | integer ) >>
                                rootNode( ch_p( 'e' ) | ch_p( 'E' ) ) >>
                                ( ch_p( '-' ) >> integer |
                                  discard_node_d[ ch_p( '+' ) ] >> integer |
                                  integer );

                negative    = rootNode( ch_p( '-' ) ) >> factor;

                identifier  =   leafNode( alpha_p >> *( alnum_p | ch_p( '_' ) ) );

                variable    =   identifier >> rootNode( ch_p( '.' ) ) >> identifier;

                entity_func = identifier >> entity_property >> rootNode( ch_p( '.' ) ) >> identifier;

                entity_property = + ( rootNode( ch_p( '.' ) ) >>
                                      leafNode( + ( alpha_p | ch_p( '_' ) ) ) >>
                                      discard_node_d[ ch_p( '(' ) ] >>
                                      discard_node_d[ ch_p( ')' ) ] );


                ///////////////////////////////////////////////////
                //                                               //
                //      This syntax is made such dirty syntax    //
                //      by the bug of Spirit                     //
                //                                               //
                ///////////////////////////////////////////////////

                //call_func = rootNode( +alpha_p ) >>

                call_func = (   rootNode( str_p( "eq" ) )
                                | rootNode( str_p( "neq" ) )
                                | rootNode( str_p( "gt" ) )
                                | rootNode( str_p( "lt" ) )
                                | rootNode( str_p( "geq" ) )
                                | rootNode( str_p( "leq" ) )
                                | rootNode( str_p( "and" ) )
                                | rootNode( str_p( "or" ) )
                                | rootNode( str_p( "xor" ) )
                                | rootNode( str_p( "not" ) )
                                | rootNode( str_p( "abs" ) )
                                | rootNode( str_p( "sqrt" ) )
                                | rootNode( str_p( "pow" ) )
                                | rootNode( str_p( "exp" ) )
                                | rootNode( str_p( "log10" ) )
                                | rootNode( str_p( "log" ) )
                                | rootNode( str_p( "floor" ) )
                                | rootNode( str_p( "ceil" ) )
                                | rootNode( str_p( "sin" ) )
                                | rootNode( str_p( "cos" ) )
                                | rootNode( str_p( "tan" ) )
                                | rootNode( str_p( "sinh" ) )
                                | rootNode( str_p( "cosh" ) )
                                | rootNode( str_p( "tanh" ) )
                                | rootNode( str_p( "asin" ) )
                                | rootNode( str_p( "acos" ) )
                                | rootNode( str_p( "atan" ) )
                                | rootNode( str_p( "fact" ) )
                                | rootNode( str_p( "asinh" ) )
                                | rootNode( str_p( "acosh" ) )
                                | rootNode( str_p( "atanh" ) )
                                | rootNode( str_p( "asech" ) )
                                | rootNode( str_p( "acsch" ) )
                                | rootNode( str_p( "acoth" ) )
                                | rootNode( str_p( "sech" ) )
                                | rootNode( str_p( "csch" ) )
                                | rootNode( str_p( "coth" ) )
                                | rootNode( str_p( "asec" ) )
                                | rootNode( str_p( "acsc" ) )
                                | rootNode( str_p( "acot" ) )
                                | rootNode( str_p( "sec" ) )
                                | rootNode( str_p( "csc" ) )
                                | rootNode( str_p( "cot" ) )
                            ) >>
                            inner_node_d[ ch_p( '(' ) >>
                                          ( expression >>
                                            *( discard_node_d[ ch_p( ',' ) ] >>
                                               expression ) ) >>
                                          ch_p( ')' ) ];


                group       =   inner_node_d[ ch_p( '(' ) >> expression >> ch_p( ')' )];

                constant    =   exponent | floating | integer;

                factor      =   call_func
                                |   entity_func
                                |   variable
                                |   constant
                                |   group
                                |   identifier
                                |   negative;

                power = factor >> *( rootNode( ch_p( '^' ) ) >> factor );

                term        =  power >>
                               *( ( rootNode( ch_p( '*' ) ) >> power )
                                  |  ( rootNode( ch_p( '/' ) ) >> power ) );
                //|  ( rootNode( ch_p('^') ) >> power ) );


                expression  =  term >>
                               *( ( rootNode( ch_p( '+' ) ) >> term )
                                  |  ( rootNode( ch_p( '-' ) ) >> term ) );
            }

            rule<ScannerT, PARSER_CONTEXT, parser_tag<VARIABLE> >     variable;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<CALL_FUNC> >    call_func;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<EXPRESSION> >   expression;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<TERM> >         term;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<POWER> >        power;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<FACTOR> >       factor;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<FLOAT> >        floating;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<EXPONENT> >     exponent;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<INTEGER> >      integer;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<NEGATIVE> >     negative;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<GROUP> >        group;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<IDENTIFIER> >   identifier;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<CONSTANT> >     constant;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<ENTITY_FUNC> >  entity_func;
            rule<ScannerT, PARSER_CONTEXT, parser_tag<ENTITY_PROPERTY> >
            entity_property;

            rule<ScannerT, PARSER_CONTEXT, parser_tag<EXPRESSION> > const&
            start() const
            {
                return expression;
            }
        };

#undef leafNode
#undef rootNode

    };

public:

    ExpressionCompiler( Process* process, PropertyMap* propertyMap )
            :
            processPtr( process ),
            propertyMapPtr( propertyMap )
    {
        if ( this->constantMap.empty() == true ||
                this->functionMap1.empty() == true )
        {
            fillMap();
        }
    }


    ~ExpressionCompiler()
    {
        ; // do nothing
    }

    typedef char const*         iterator_t;
    typedef tree_match<iterator_t> parse_tree_match_t;
    typedef parse_tree_match_t::tree_iterator TreeIterator;

    const Code compileExpression( const String& expression );

protected:

    template < class INSTRUCTION >
    static void appendInstruction( Code& code,
                                   const INSTRUCTION& instruction )
    {
        Code::size_type codeSize( code.size() );
        code.resize( codeSize + sizeof( INSTRUCTION ) );
        new ( &code[codeSize] ) INSTRUCTION( instruction );
    }

    /**template < Opcode OPCODE >
       static void appendSimpleInstruction( Code& code )
       {
       appendInstruction( code, Instruction<OPCODE>() );
       }*/

    static void
    appendVariableReferenceMethodInstruction( Code& code,
            VariableReference*              variableReference,
            const String& methodName );

    static void
    appendEntityMethodInstruction( Code& code,
                                   Entity* entityPtr,
                                   const String& methodName );


private:

    static void fillMap();

    void compileTree( TreeIterator const& treeIterator, Code& code );
    void compileEntityProperty
    ( TreeIterator const& treeIterator, Code& code,
      Entity* entityPtr, const String methodName );

    void throw_exception( String type, String string );

private:

    Process*      processPtr;
    PropertyMap*  propertyMapPtr;

    static ConstantMap        constantMap;
    static FunctionMap1       functionMap1;
    static FunctionMap2       functionMap2;

}
; // ExpressionCompiler





template <>
class ExpressionCompiler::InstructionBase<ExpressionCompiler::NoOperand>
            :
            public ExpressionCompiler::InstructionHead
{
public:

    InstructionBase( Opcode opcode )
            :
            InstructionHead( opcode )
    {
        ; // do nothing
    }

    InstructionBase( Opcode, const NoOperand& );

};




#define SPECIALIZE_OPCODE2OPERAND( OP, OPE )\
  template<> class ExpressionCompiler::Opcode2Operand<ExpressionCompiler::OP>\
  {\
  public:\
    typedef ExpressionCompiler::OPE type;\
  };



SPECIALIZE_OPCODE2OPERAND( PUSH_REAL,                Real );
//SPECIALIZE_OPCODE2OPERAND( PUSH_POINTER,             Pointer );
SPECIALIZE_OPCODE2OPERAND( LOAD_REAL,                Real* const );
//SPECIALIZE_OPCODE2OPERAND( CALL_FUNC0,           RealFunc0 );
SPECIALIZE_OPCODE2OPERAND( CALL_FUNC1,               RealFunc1 );
SPECIALIZE_OPCODE2OPERAND( CALL_FUNC2,               RealFunc2 );
SPECIALIZE_OPCODE2OPERAND( OBJECT_METHOD_REAL,
                           RealObjectMethodProxy );
SPECIALIZE_OPCODE2OPERAND( OBJECT_METHOD_INTEGER,
                           IntegerObjectMethodProxy );


#define DEFINE_OPCODE2INSTRUCTION( CODE )\
  template<> class\
    ExpressionCompiler::Opcode2Instruction<ExpressionCompiler::CODE>\
  {\
  public:\
    typedef ExpressionCompiler::Instruction<ExpressionCompiler::CODE> type;\
    typedef type::Operand operandtype;\
  }


DEFINE_OPCODE2INSTRUCTION( PUSH_REAL );
//DEFINE_OPCODE2INSTRUCTION( PUSH_POINTER );
DEFINE_OPCODE2INSTRUCTION( NEG );
DEFINE_OPCODE2INSTRUCTION( ADD );
DEFINE_OPCODE2INSTRUCTION( SUB );
DEFINE_OPCODE2INSTRUCTION( MUL );
DEFINE_OPCODE2INSTRUCTION( DIV );
//DEFINE_OPCODE2INSTRUCTION( POW );
DEFINE_OPCODE2INSTRUCTION( LOAD_REAL );
//DEFINE_OPCODE2INSTRUCTION( CALL_FUNC0 );
DEFINE_OPCODE2INSTRUCTION( CALL_FUNC1 );
DEFINE_OPCODE2INSTRUCTION( CALL_FUNC2 );
DEFINE_OPCODE2INSTRUCTION( OBJECT_METHOD_INTEGER );
DEFINE_OPCODE2INSTRUCTION( OBJECT_METHOD_REAL );
DEFINE_OPCODE2INSTRUCTION( RET );


const ExpressionCompiler::Code
ExpressionCompiler::compileExpression( const String& expression )
{
    Code code;
    CompileGrammar grammer;

    tree_parse_info<>
    info( ast_parse( expression.c_str(), grammer, space_p ) );

    if ( expression.length() == 0 )
    {
        THROW_EXCEPTION( UnexpectedError,
                         "Expression is empty\nClass : " +
                         String( this->processPtr->getClassName() ) +
                         "\nProcessID : " + String( this->processPtr->getID() ) );
    }

    else
    {
        if ( info.full )
        {
            compileTree( info.trees.begin(), code );

            // place RET at the tail.
            appendInstruction( code, Instruction<RET>() );
        }
        else
        {
            THROW_EXCEPTION( UnexpectedError,
                             "Parse error in the expression.\nExpression : "
                             + expression + "\nClass : "
                             + String( this->processPtr->getClassName() )
                             + "\nProcessID : "
                             + String( this->processPtr->getID() ) );
        }
    }

    return code;
}


void ExpressionCompiler::fillMap()
{

    // set ConstantMap
    constantMap["true"]  = 1.0;
    constantMap["false"] = 0.0;
    constantMap["pi"]    = M_PI;
    constantMap["Nn"]   = std::numeric_limits<Real>::quiet_NaN();
    constantMap["INF"]   = std::numeric_limits<Real>::infinity();
    constantMap["N_A"]   = N_A;
    constantMap["exp"]   = M_E;


    // set FunctionMap1
    functionMap1["abs"]   = std::fabs;
    functionMap1["sqrt"]  = std::sqrt;
    functionMap1["exp"]   = std::exp;
    functionMap1["log10"] = std::log10;
    functionMap1["log"]   = std::log;
    functionMap1["floor"] = std::floor;
    functionMap1["ceil"]  = std::ceil;
    functionMap1["sin"]   = std::sin;
    functionMap1["cos"]   = std::cos;
    functionMap1["tan"]   = std::tan;
    functionMap1["sinh"]  = std::sinh;
    functionMap1["cosh"]  = std::cosh;
    functionMap1["tanh"]  = std::tanh;
    functionMap1["asin"]  = std::asin;
    functionMap1["acos"]  = std::acos;
    functionMap1["atan"]  = std::atan;
    functionMap1["fact"]  = fact;
    functionMap1["asinh"] = asinh;
    functionMap1["acosh"] = acosh;
    functionMap1["atanh"] = atanh;
    functionMap1["asech"] = asech;
    functionMap1["acsch"] = acsch;
    functionMap1["acoth"] = acoth;
    functionMap1["sech"]  = sech;
    functionMap1["csch"]  = csch;
    functionMap1["coth"]  = coth;
    functionMap1["asec"]  = asec;
    functionMap1["acsc"]  = acsc;
    functionMap1["acot"]  = acot;
    functionMap1["sec"]   = sec;
    functionMap1["csc"]   = csc;
    functionMap1["cot"]   = cot;
    functionMap1["not"]   = libecs::real_not;


    // set FunctionMap2
    functionMap2["pow"]   = pow;
    functionMap2["and"]   = libecs::real_and;
    functionMap2["or"]    = libecs::real_or;
    functionMap2["xor"]   = libecs::real_xor;
    functionMap2["eq"]    = libecs::real_eq;
    functionMap2["neq"]   = libecs::real_neq;
    functionMap2["gt"]    = libecs::real_gt;
    functionMap2["lt"]    = libecs::real_lt;
    functionMap2["geq"]   = libecs::real_geq;
    functionMap2["leq"]   = libecs::real_leq;

}


#define APPEND_OBJECT_METHOD_REAL( OBJECT, CLASSNAME, METHODNAME )\
 appendInstruction\
   ( code, \
     Instruction<OBJECT_METHOD_REAL>\
     ( RealObjectMethodProxy::\
       create< CLASSNAME, & CLASSNAME::METHODNAME >\
       ( OBJECT ) ) ) // \

#define APPEND_OBJECT_METHOD_INTEGER( OBJECT, CLASSNAME, METHODNAME )\
 appendInstruction\
   ( code, \
     Instruction<OBJECT_METHOD_INTEGER>\
     ( IntegerObjectMethodProxy::\
       create< CLASSNAME, & CLASSNAME::METHODNAME >\
       ( OBJECT ) ) ) // \


void
ExpressionCompiler::
appendVariableReferenceMethodInstruction( Code& code,
        VariableReference*          variableReference,
        const String& methodName )
{

    if ( methodName == "Value" )
    {
        APPEND_OBJECT_METHOD_REAL( variableReference, VariableReference,
                                   getValue );
    }
    else if ( methodName == "Velocity" )
    {
        APPEND_OBJECT_METHOD_REAL( variableReference, VariableReference,
                                   getVelocity );
    }
    else if ( methodName == "Coefficient" )
    {
        APPEND_OBJECT_METHOD_INTEGER( variableReference, VariableReference,
                                      getCoefficient );
    }

    /**else if( str_child2 == "Fixed" ){
       code.push_back(
       new OBJECT_METHOD_REAL( variableReference,
       &libecs::VariableReference::isFixed ) );
       }*/

    else
    {
        THROW_EXCEPTION
        ( NotFound,
          "VariableReference attribute [" +
          methodName + "] not found." );
    }


}

void
ExpressionCompiler::
appendEntityMethodInstruction( Code& code,
                               Entity* entityPtr,
                               const String& methodName )
{
    
    THROW_EXCEPTION
        ( NotFound,
          "Entity attribute [" +
          methodName + "] not found." );


}

#undef APPEND_OBJECT_METHOD_REAL
#undef APPEND_OBJECT_METHOD_INTEGER



void
ExpressionCompiler::throw_exception( String exceptionType,
                                     String exceptionString )
{
    if ( exceptionType == "UnexpeptedError" )
    {
        THROW_EXCEPTION( UnexpectedError, exceptionString );
    }
    else if ( exceptionType == "NoSlot" )
    {
        THROW_EXCEPTION( NoSlot, exceptionString );
    }
    else if ( exceptionType == "NotFound" )
    {
        THROW_EXCEPTION( NotFound, exceptionString );
    }
    else
    {
        THROW_EXCEPTION( UnexpectedError, exceptionString );
    }
}

/**
   This function is ExpressionCompiler subclass member function.
   This member function evaluates AST tree and makes binary codes.
*/

void ExpressionCompiler::compileTree( TreeIterator const& treeIterator,
                                      Code& code )
{
    /**
       compile AST
    */

    switch ( treeIterator->value.id().to_long() )
    {
        /**
        Floating Grammar compile
        */

    case CompileGrammar::FLOAT :
    {
        assert( treeIterator->children.size() == 0 );

        const String floatString( treeIterator->value.begin(),
                                  treeIterator->value.end() );

        const Real floatValue = stringCast<Real>( floatString );

        appendInstruction( code, Instruction<PUSH_REAL>( floatValue ) );

        return;
    }

    /**
       Integer Grammar compile
    */

    case CompileGrammar::INTEGER :
    {
        assert( treeIterator->children.size() == 0 );

        const String integerString( treeIterator->value.begin(),
                                    treeIterator->value.end() );

        const Real integerValue = stringCast<Real>( integerString );

        appendInstruction( code, Instruction<PUSH_REAL>( integerValue ) );

        return;

    }

    /**
       Grammar compile
    */

    case CompileGrammar::EXPONENT:
    {
        assert( *treeIterator->value.begin() == 'E' ||
                *treeIterator->value.begin() == 'e' );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );

        const String baseString( childTreeIterator->value.begin(),
                                 childTreeIterator->value.end() );

        const String
        exponentString( ( childTreeIterator + 1 )->value.begin(),
                        ( childTreeIterator + 1 )->value.end() );

        const Real baseValue = stringCast<Real>( baseString );

        if ( exponentString != "-" )
        {
            const Real
            exponentValue = stringCast<Real>( exponentString );

            appendInstruction
            ( code, Instruction<PUSH_REAL>
              ( baseValue * pow( 10, exponentValue ) ) );
        }
        else
        {
            const String
            exponentString1( ( childTreeIterator + 2 )->value.begin(),
                             ( childTreeIterator + 2 )->value.end() );

            const Real
            exponentValue = stringCast<Real>( exponentString1 );

            appendInstruction
            ( code,
              Instruction<PUSH_REAL>
              ( baseValue * pow( 10, -exponentValue ) ) );
        }

        return;
    }



    /**
    Call_Func Grammar compile
    */

    case CompileGrammar::CALL_FUNC :
    {
        Integer childTreeSize( treeIterator->children.size() );

        const String functionString( treeIterator->value.begin(),
                                     treeIterator->value.end() );


        assert( childTreeSize != 0 );

        FunctionMap1::iterator functionMap1Iterator;
        FunctionMap2::iterator functionMap2Iterator;


        if ( childTreeSize == 1 )
        {
            functionMap1Iterator =
                this->functionMap1.find( functionString );

            TreeIterator const&
            childTreeIterator( treeIterator->children.begin() );


            if ( childTreeIterator->value.id() == CompileGrammar::INTEGER ||
                    childTreeIterator->value.id() == CompileGrammar::FLOAT )
            {
                const String
                argumentString( childTreeIterator->value.begin(),
                                childTreeIterator->value.end() );

                const Real
                argumentValue = stringCast<Real>( argumentString );

                if ( functionMap1Iterator != this->functionMap1.end() )
                {
                    appendInstruction
                    ( code, Instruction<PUSH_REAL>
                      ( ( *functionMap1Iterator->second )
                        ( argumentValue ) ) );
                }
                else
                {
                    functionMap2Iterator =
                        this->functionMap2.find( functionString );

                    if ( functionMap2Iterator != this->functionMap2.end() )
                    {
                        ExpressionCompiler::throw_exception
                        ( "UnexpectedError",
                          "[ " + functionString +
                          " ] function. Too few arguments\nProcessID : "
                          + this->processPtr->getID() );
                    }
                    else
                    {
                        ExpressionCompiler::throw_exception
                        ( "NoSlot",
                          "[ " + functionString +
                          String( " ] : No such function." ) +
                          "\nProcessID : " + this->processPtr->getID() );
                    }
                }
            }
            else
            {
                compileTree( childTreeIterator, code );

                if ( functionMap1Iterator != this->functionMap1.end() )
                {
                    appendInstruction
                    ( code, Instruction<CALL_FUNC1>
                      ( functionMap1Iterator->second ) );
                }
                else
                {
                    functionMap2Iterator =
                        this->functionMap2.find( functionString );

                    if ( functionMap2Iterator != this->functionMap2.end() )
                    {
                        ExpressionCompiler::throw_exception
                        ( "UnexpectedError",
                          "[ " + functionString +
                          " ] function. Too few arguments\nProcessID : "
                          + this->processPtr->getID() );
                    }
                    else
                    {
                        ExpressionCompiler::throw_exception
                        ( "NoSlot",
                          "[ " + functionString +
                          String( " ] : No such function." ) +
                          "\nProcessID : " + this->processPtr->getID() );
                    }
                }
            }
        }

        else if ( childTreeSize == 2 )
        {
            TreeIterator const&
            childTreeIterator( treeIterator->children.begin() );

            compileTree( childTreeIterator, code );
            compileTree( childTreeIterator + 1, code );


            functionMap2Iterator =
                this->functionMap2.find( functionString );

            if ( functionMap2Iterator != this->functionMap2.end() )
            {
                appendInstruction
                ( code, Instruction<CALL_FUNC2>
                  ( functionMap2Iterator->second ) );
            }
            else
            {
                functionMap1Iterator =
                    this->functionMap1.find( functionString );

                if ( functionMap1Iterator != this->functionMap1.end() )
                {
                    ExpressionCompiler::throw_exception
                    ( "UnexpectedError",
                      "[ " + functionString +
                      " ] function. Too many arguments\nProcessID : " +
                      this->processPtr->getID() );
                }
                else
                {
                    ExpressionCompiler::throw_exception
                    ( "NotFound",
                      "[ " + functionString +
                      String( " ] : No such function." ) +
                      "\nProcessID : " +
                      this->processPtr->getID() );
                }
            }
        }

        else
        {
            ExpressionCompiler::throw_exception
            ( "UnexpectedError",
              " : Too many arguments\nProcessID : " +
              this->processPtr->getID() );
        }

        return;
    }


    /**
       Entity_Func Grammar compile
    */

    case CompileGrammar::ENTITY_FUNC :
    {
        assert( treeIterator->children.size() >= 3 );
        Integer childTreeSize( treeIterator->children.size() );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );

        const String classString( childTreeIterator->value.begin(),
                                  childTreeIterator->value.end() );

        assert( *treeIterator->value.begin() == '.' );

        if ( classString == "self" ) // Process Class
        {
            Entity* entityPtr( this->processPtr->getSuperEntity() );

            const String methodName
            ( ( childTreeIterator + childTreeSize - 1 )->value.begin(),
              ( childTreeIterator + childTreeSize - 1 )->value.end() );

            compileEntityProperty( childTreeIterator + 1,
                                   code,
                                   entityPtr,
                                   methodName );
        }

        else // VariableReference Class
        {
            const VariableReference&              variableReference( this->processPtr->
                    getVariableReference( classString ) );

            Entity* const entityPtr( variableReference.getSuperEntity() );

            const String methodName
            ( ( childTreeIterator + childTreeSize - 1 )->value.begin(),
              ( childTreeIterator + childTreeSize - 1 )->value.end() );

            compileEntityProperty( childTreeIterator + 1,
                                   code,
                                   entityPtr,
                                   methodName );
        }
        return;
    }


    /**
       Variable Grammar compile
    */

    case CompileGrammar::VARIABLE :
    {
        assert( *treeIterator->value.begin() == '.' );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );

        const String
        variableReferenceString( childTreeIterator->value.begin(),
                                 childTreeIterator->value.end() );

        const String
        variableReferenceMethodString
        ( ( childTreeIterator + 1 )->value.begin(),
          ( childTreeIterator + 1 )->value.end() );

        const VariableReference&          variableReference
        ( this->processPtr->
          getVariableReference( variableReferenceString ) );

        appendVariableReferenceMethodInstruction
        ( code,
          const_cast<VariableReference*>( &variableReference ),
          variableReferenceMethodString );

        return;

    }



    /**
       Identifier Grammar compile
    */

    case CompileGrammar::IDENTIFIER :
    {
        assert( treeIterator->children.size() == 0 );

        const String identifierString( treeIterator->value.begin(),
                                       treeIterator->value.end() );

        ConstantMap::iterator constantMapIterator;
        PropertyMap::iterator propertyMapIterator;

        constantMapIterator =
            this->constantMap.find( identifierString );
        propertyMapIterator =
            this->propertyMapPtr->find( identifierString );


        if ( constantMapIterator != this->constantMap.end() )
        {
            appendInstruction
            ( code,
              Instruction<PUSH_REAL>( constantMapIterator->second ) );
        }

        else if ( propertyMapIterator != this->propertyMapPtr->end() )
        {
            appendInstruction
            ( code, Instruction<LOAD_REAL>
              ( &( propertyMapIterator->second ) ) );
        }

        else
        {
            ExpressionCompiler::throw_exception
            ( "NotFound",
              "[ " + identifierString +
              " ] No such Property slot.\nProcessID : "
              + this->processPtr->getID() );
        }

        return;
    }



    /**
       Negative Grammar compile 
    */

    case CompileGrammar::NEGATIVE :
    {
        assert( *treeIterator->value.begin() == '-' );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );


        if ( childTreeIterator->value.id() == CompileGrammar::INTEGER ||
                childTreeIterator->value.id() == CompileGrammar::FLOAT )
        {
            const String
            valueString( childTreeIterator->value.begin(),
                         childTreeIterator->value.end() );

            const Real
            value = stringCast<Real>( valueString );

            appendInstruction( code, Instruction<PUSH_REAL>( -value ) );
        }
        else
        {
            compileTree( childTreeIterator, code );

            appendInstruction( code, Instruction<NEG>() );
        }

        return;

    }



    /**
       Power Grammar compile
    */

    case CompileGrammar::POWER :
    {
        assert( treeIterator->children.size() == 2 );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );


        if ( ( childTreeIterator->value.id() == CompileGrammar::INTEGER ||
                childTreeIterator->value.id() == CompileGrammar::FLOAT ) &&
                ( ( childTreeIterator + 1 )->value.id() == CompileGrammar::INTEGER ||
                  ( childTreeIterator + 1 )->value.id() == CompileGrammar::FLOAT ) )
        {

            const String
            argumentString1( childTreeIterator->value.begin(),
                             childTreeIterator->value.end() );

            const String
            argumentString2( ( childTreeIterator + 1 )->value.begin(),
                             ( childTreeIterator + 1 )->value.end() );

            const Real
            argumentValue1 = stringCast<Real>( argumentString1 );
            const Real
            argumentValue2 = stringCast<Real>( argumentString2 );


            if ( *treeIterator->value.begin() == '^' )
            {
                appendInstruction
                ( code, Instruction<PUSH_REAL>
                  ( pow( argumentValue1, argumentValue2 ) ) );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalid operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }

            return;
        }
        else
        {
            compileTree( treeIterator->children.begin(), code );
            compileTree( treeIterator->children.begin() + 1, code );

            if ( *treeIterator->value.begin() == '^' )
            {
                RealFunc2 powFunc( this->functionMap2.find( "pow" )->second );
                appendInstruction( code,
                                   Instruction<CALL_FUNC2>( powFunc ) );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalud operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }

            return;
        }

        return;

    }



    /**
       Term Grammar compile
    */

    case CompileGrammar::TERM :
    {

        assert( treeIterator->children.size() == 2 );


        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );


        if ( ( childTreeIterator->value.id() == CompileGrammar::INTEGER ||
                childTreeIterator->value.id() == CompileGrammar::FLOAT ) &&
                ( ( childTreeIterator + 1 )->value.id() == CompileGrammar::INTEGER ||
                  ( childTreeIterator + 1 )->value.id() == CompileGrammar::FLOAT ) )
        {

            const String term1String( childTreeIterator->value.begin(),
                                      childTreeIterator->value.end() );

            const String
            term2String( ( childTreeIterator + 1 )->value.begin(),
                         ( childTreeIterator + 1 )->value.end() );

            const Real term1Value = stringCast<Real>( term1String );
            const Real term2Value = stringCast<Real>( term2String );


            if ( *treeIterator->value.begin() == '*' )
            {
                appendInstruction
                ( code,
                  Instruction<PUSH_REAL>( term1Value * term2Value ) );
            }

            else if ( *treeIterator->value.begin() == '/' )
            {
                appendInstruction
                ( code,
                  Instruction<PUSH_REAL>( term1Value / term2Value ) );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalid operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }

            return;
        }
        else
        {
            compileTree( childTreeIterator, code );
            compileTree( childTreeIterator + 1, code );

            if ( *treeIterator->value.begin() == '*' )
            {
                appendInstruction( code, Instruction<MUL>() );
            }

            else if ( *treeIterator->value.begin() == '/' )
            {
                appendInstruction( code, Instruction<DIV>() );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalid operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }

            return;
        }

        return;

    }



    /**
       Expression Grammar compile
    */

    case CompileGrammar::EXPRESSION :
    {

        assert( treeIterator->children.size() == 2 );

        TreeIterator const&
        childTreeIterator( treeIterator->children.begin() );


        if ( ( childTreeIterator->value.id() == CompileGrammar::INTEGER ||
                childTreeIterator->value.id() == CompileGrammar::FLOAT ) &&
                ( ( childTreeIterator + 1 )->value.id() == CompileGrammar::INTEGER ||
                  ( childTreeIterator + 1 )->value.id() == CompileGrammar::FLOAT ) )
        {
            const String term1String( childTreeIterator->value.begin(),
                                      childTreeIterator->value.end() );

            const String
            term2String( ( childTreeIterator + 1 )->value.begin(),
                         ( childTreeIterator + 1 )->value.end() );

            const Real term1Value = stringCast<Real>( term1String );
            const Real term2Value = stringCast<Real>( term2String );


            if ( *treeIterator->value.begin() == '+' )
            {
                appendInstruction
                ( code,
                  Instruction<PUSH_REAL>( term1Value + term2Value ) );
            }

            else if ( *treeIterator->value.begin() == '-' )
            {
                appendInstruction
                ( code,
                  Instruction<PUSH_REAL>( term1Value - term2Value ) );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalid operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }
        }
        else
        {
            compileTree( childTreeIterator, code );
            compileTree( childTreeIterator + 1, code );


            if ( *treeIterator->value.begin() == '+' )
            {
                appendInstruction( code, Instruction<ADD>() );
            }

            else if ( *treeIterator->value.begin() == '-' )
            {
                appendInstruction( code, Instruction<SUB>() );
            }

            else
            {
                ExpressionCompiler::throw_exception
                ( "UnexpectedError",
                  String( "Invalid operation" ) +
                  "\nProcessID : " + this->processPtr->getID() );
            }
        }

        return;

    }

    default :
    {
        ExpressionCompiler::throw_exception
        ( "UnexpectedError",
          "syntax error.\nProcessID : " + this->processPtr->getID() );

        return;
    }
    }
}


void
ExpressionCompiler::compileEntityProperty( TreeIterator const& treeIterator,
        Code& code,
        Entity* entityPtr,
        const String methodName )
{
    TreeIterator const&
    childTreeIterator( treeIterator->children.begin() );

    const String childString( childTreeIterator->value.begin(),
                              childTreeIterator->value.end() );

    assert( *treeIterator->value.begin() == '.' );

    if ( childString == "getSuperEntity" )
    {
        appendEntityMethodInstruction( code,
                                       entityPtr,
                                       methodName );
    }
    else if ( childString == "." )
    {
        Entity* entityPtr( entityPtr->getSuperEntity() );

        compileEntityProperty( childTreeIterator,
                               code,
                               entityPtr,
                               methodName );
    }
    else
    {
        ExpressionCompiler::throw_exception
        ( "UnexpectedError",
          String( "Entity function parse error" ) +
          "\nProcessID : " + this->processPtr->getID() );
    }
}



// this should be moved to .cpp

ConstantMap ExpressionCompiler::constantMap;
FunctionMap1 ExpressionCompiler::functionMap1;
FunctionMap2 ExpressionCompiler::functionMap2;

#endif /* __EXPRESSIONCOMPILER_HPP */


