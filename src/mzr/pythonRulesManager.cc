
#include "pythonRulesManager.hh"


#include <iostream>

namespace mzr
{

  void PythonRulesManager::DEBUG_doInterestingStuff()
  {
    
    PyObject_CallMethod(mzrFileConverterClassInst, "DEBUGPRINT", NULL);
    return;
  }

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
    strcpy( explicitSpeciesPythonFunctionName, "addExplicitSpeciesStatement" );
    strcpy (singleStringTuple, "s");

    PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pFileName;

    Py_Initialize();

    pName = PyString_FromString( "moleculizer" );
    pFileName = PyString_FromString( "dmp_tmp_out" );

    pModule = PyImport_Import( pName );
    pDict = PyModule_GetDict( pModule );
    pFunc = PyDict_GetItemString( pDict, "MoleculizerRulesFile" );

    pArgs = PyTuple_New(1);
    PyTuple_SetItem(pArgs, 0, pFileName);

    mzrFileConverterClassInst = PyObject_CallObject(pFunc, pArgs);

    this->DEBUG_doInterestingStuff();
  }

  std::string
  PythonRulesManager::getXmlString() const
  {

    PyErr_Print();

    PyObject* fileAsString;

    if (!isInitialized() ) throw std::exception();
    fileAsString = PyObject_CallMethod(mzrFileConverterClassInst, getFileStringFunctionName, NULL);

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
