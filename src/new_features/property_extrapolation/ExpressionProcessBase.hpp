//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//       This file is part of the E-Cell System
//
//       Copyright (C) 1996-2008 Keio University
//       Copyright (C) 2005-2008 The Molecular Sciences Institute
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell System is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// E-Cell System is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with E-Cell System -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//END_HEADER
//
// authors:
//   Koichi Takahashi
//   Tatsuya Ishida
//
// E-Cell Project.
//

#ifndef __EXPRESSIONPROCESSBASE_HPP
#define __EXPRESSIONPROCESSBASE_HPP

#define EXPRESSION_PROCESS_USE_JIT 0

#define ENABLE_STACKOPS_FOLDING 1


#include <cassert>
#include <limits>

#include "ExpressionCompiler.hpp"

//#if defined( EXPRESSIONPROCESS_USE_JIT )
//#include "JITExpressionProcessBase"
//#else /* defined( EXPRESSIONPROCESS_USE_JIT ) */
//#include "SVMExpressionProcessBase"
//#endif /* defined( EXPRESSIONPROCESS_USE_JIT ) */


USE_LIBECS;

LIBECS_DM_CLASS( ExpressionProcessBase, Process )
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
        {
            // ; do nothing
        }

        ~VirtualMachine()
        {}

        const Real execute( const Code& code );

    };


public:

    LIBECS_DM_OBJECT_ABSTRACT( ExpressionProcessBase )
    {
        INHERIT_PROPERTIES( Process );

        PROPERTYSLOT_SET_GET( String, Expression );
    }


    ExpressionProcessBase()
            :
            recompileFlag( true )
    {
        // ; do nothing
    }

    virtual ~ExpressionProcessBase()
    {
        // ; do nothing
    }

    SET_METHOD( String, Expression )
    {
        this->expression = value;
        this->recompileFlag = true;
    }

    GET_METHOD( String, Expression )
    {
        return this->expression;
    }


    void defaultSetProperty( const String& propertyName,
                             const Polymorph& value )
    {
        this->propertyMap[ propertyName ] = value.asReal();
    }

    const Polymorph defaultGetProperty( const String& propertyName ) const
    {
        PropertyMap::const_iterator          propertyMapIterator( this->propertyMap.find( propertyName ) );

        if ( propertyMapIterator != this->propertyMap.end() )
        {
            return propertyMapIterator->second;
        }
        else
        {
            THROW_EXCEPTION( NoSlot, getClassNameString() +
                             " : Property [" + propertyName +
                             "] is not defined " );
        }
    }

    const Polymorph defaultGetPropertyList() const
    {
        PolymorphVector vector;

        for ( PropertyMap::const_iterator                  propertyMapIterator( this->propertyMap.begin() );
                propertyMapIterator != this->propertyMap.end();
                ++propertyMapIterator )
        {
            vector.push_back( propertyMapIterator->first );
        }

        return vector;
    }

    const Polymorph
    defaultGetPropertyAttributes( const String& propertyName ) const
    {
        PolymorphVector vector;

        Integer propertyFlag( 1 );

        vector.push_back( propertyFlag ); // isSetable
        vector.push_back( propertyFlag ); // isGetable
        vector.push_back( propertyFlag ); // isLoadable
        vector.push_back( propertyFlag ); // isSavable

        return vector;
    }


    void compileExpression()
    {
        ExpressionCompiler compiler( this, &( getPropertyMap() ) );

        this->compiledCode.clear();
        this->compiledCode = compiler.compileExpression( this->expression );

        // virtualMachine.resize( compiler.getStackSize() );
    }

    const PropertyMap& getPropertyMap() const
    {
        return this->propertyMap;
    }

    virtual void initialize()
    {
        Process::initialize();

        if ( this->recompileFlag )
        {
            compileExpression();
            this->recompileFlag = false;
        }
    }

protected:

    PropertyMap& getPropertyMap()
    {
        return this->propertyMap;
    }


protected:

    String    expression;

    Code compiledCode;
    VirtualMachine virtualMachine;

    bool recompileFlag;

    PropertyMap propertyMap;
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
            = ( instruction->getOperand() )( ( stackPtr - 1 )->real,
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

            bypass = static_cast<Real>( ( instruction->getOperand() )() );

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


LIBECS_DM_INIT_STATIC( ExpressionProcessBase, Process );

#endif /* __EXPRESSIONPROCESSBASE_HPP */


