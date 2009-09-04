//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of Libmoleculizer
//
//        Copyright (C) 2001-2009 The Molecular Sciences Institute.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Moleculizer is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// Moleculizer is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Moleculizer; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307,  USA
//
// END HEADER
//
// Original Author:
//   Nathan Addy, Scientific Programmer, Molecular Sciences Institute, 2009
//
// Modifing Authors:
//
//

#include "pythonRulesManager.hh"
#include "mzrException.hh"

namespace mzr
{

  PythonRulesManager::~PythonRulesManager()
  {
    delete paramPythonFunctionName;
    delete molPythonFunctionName;
    delete alloPlexPythonFunctionName;
    delete alloOmniPythonFunctionName;
    delete dimerGenPythonFunctionName;
    delete omniGenPythonFunctionName;
    delete uniMolGenPythonFunctionName;
    delete speciesStreamPythonFunctionName;
    delete singleStringTuple;
    delete getFileStringFunctionName;
    delete addWholeRulesFileFunctionName;
    delete addWholeRulesStringFunctionName;
    
    // PyDecref( mzrFileConverterClassInst );
  }


  PythonRulesManager::PythonRulesManager()
    :
    _isInitialized( true )
  {
    paramPythonFunctionName = new char[256];
    modPythonFunctionName = new char[256];
    molPythonFunctionName = new char[256];
    alloPlexPythonFunctionName = new char[256];
    alloOmniPythonFunctionName = new char[256];
    dimerGenPythonFunctionName = new char[256];
    omniGenPythonFunctionName = new char[256];
    uniMolGenPythonFunctionName = new char[256];
    speciesStreamPythonFunctionName = new char[256];
    singleStringTuple = new char[256];
    getFileStringFunctionName = new char[256];
    addWholeRulesFileFunctionName = new char[256];
    addWholeRulesStringFunctionName = new char[256];
    explicitSpeciesPythonFunctionName = new char[256];

    strcpy( paramPythonFunctionName, "addParameterStatement");
    strcpy( modPythonFunctionName, "addModificationStatement");
    strcpy( molPythonFunctionName, "addMolStatement");
    strcpy( alloPlexPythonFunctionName, "addAllostericPlexStatement");
    strcpy( alloOmniPythonFunctionName, "addOmniGenStatement");
    strcpy( dimerGenPythonFunctionName, "addDimerizationGenStatement");
    strcpy( omniGenPythonFunctionName, "addOmniGenStatement");
    strcpy( uniMolGenPythonFunctionName, "addUniMolGenStatement");
    strcpy( speciesStreamPythonFunctionName, "addSpeciesStreamStatement");
    strcpy( getFileStringFunctionName, "writeToString");
    strcpy( addWholeRulesFileFunctionName, "addWholeRulesFile");
    strcpy( addWholeRulesStringFunctionName, "addWholeRulesString");
    strcpy( explicitSpeciesPythonFunctionName, "addExplicitSpeciesStatement" );
    strcpy (singleStringTuple, "s");

    PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pFileName;

    Py_Initialize();

    pName = PyString_FromString( "moleculizer" );

    pModule = PyImport_Import( pName );

    if (!pModule)
    {
	_isInitialized = false;
	PyErr_Print();
	throw mzrPythonXcpt("The moleculizer module could not be loaded.  Please check to make sure it is installed or contact your system administrator.");
    }

    pDict = PyModule_GetDict( pModule );
    pFunc = PyDict_GetItemString( pDict, "MoleculizerRulesFile" );

    if (!pModule)
    {
	_isInitialized = false;
	PyErr_Print();
	throw mzrPythonXcpt("The MoleculizerRulesFile could not be found in the moleculizer module. Please contact your system administrator");
    }

    mzrFileConverterClassInst = PyObject_CallObject(pFunc, NULL);

  }

  std::string
  PythonRulesManager::getXmlString() const
  {

    PyErr_Print();

    PyObject* fileAsString;

    if (!isInitialized() ) throw std::exception();
    fileAsString = PyObject_CallMethod(mzrFileConverterClassInst, getFileStringFunctionName, NULL);

    PyErr_Print();

    return std::string(PyString_AsString(fileAsString));
  }



  std::string 
  PythonRulesManager::addRulesFile( const std::string& rulesFile)
  {
    if( !isInitialized() ) throw std::exception();
    PyErr_Print();
    PyObject_CallMethod( mzrFileConverterClassInst, addWholeRulesFileFunctionName, singleStringTuple, rulesFile.c_str());
    PyErr_Print();
    return getXmlString();
  }


  std::string
  PythonRulesManager::addRulesString( const std::string& rulesString)
  {
    if( !isInitialized() ) throw std::exception();

    PyErr_Print();
    PyObject_CallMethod( mzrFileConverterClassInst, addWholeRulesStringFunctionName, singleStringTuple, rulesString.c_str());
    PyErr_Print();

    return getXmlString();
  }


  void PythonRulesManager::addParameterStatement(const std::string& paramLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, paramPythonFunctionName, singleStringTuple, paramLine.c_str());
  }

  void PythonRulesManager::addModificationStatement(const std::string& modLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, modPythonFunctionName, singleStringTuple, modLine.c_str());

  }

  void PythonRulesManager::addMolsStatement(const std::string& molsLIne)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, molPythonFunctionName, singleStringTuple, molsLIne.c_str());

  }

  void PythonRulesManager::addAllostericPlexStatement(const std::string& alloPlexLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, alloPlexPythonFunctionName, singleStringTuple, alloPlexLine.c_str());

  }

  void PythonRulesManager::addAllostericOmniStatement(const std::string& alloOmniLine)
  {

    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, alloOmniPythonFunctionName, singleStringTuple, alloOmniLine.c_str());

  }

  void PythonRulesManager::addDimerizationGenStatement(const std::string& dimerGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, dimerGenPythonFunctionName, singleStringTuple, dimerGenLine.c_str());
  }

  void PythonRulesManager::addOmniGenStatement(const std::string& omniGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, omniGenPythonFunctionName, singleStringTuple, omniGenLine.c_str());
  }

  void PythonRulesManager::addUniMolGenStatement(const std::string& uniMolGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, uniMolGenPythonFunctionName, singleStringTuple, uniMolGenLine.c_str());
  }

  void PythonRulesManager::addSpeciesStreamStatement(const std::string& speciesStreamLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, speciesStreamPythonFunctionName, singleStringTuple, speciesStreamLine.c_str() );
  }


}
