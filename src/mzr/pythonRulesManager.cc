
#include "pythonRulesManager.hh"


#include <iostream>

namespace mzr
{

  void PythonRulesManager::DEBUG_doInterestingStuff()
  {
    PyObject *pName, *pModule, *pDict, *pFunc, *pValue;
    PyObject *pFileName, *pArgs;

    Py_Initialize();

    pName = PyString_FromString( "moleculizer" );
    pModule = PyImport_Import( pName );
    pDict = PyModule_GetDict( pModule );
    pFunc = PyDict_GetItemString( pDict, "MoleculizerRulesFile" );
    
    pFileName = PyString_FromString("FOO");

     if (PyCallable_Check( pFunc ) )
       {

	 PyErr_Print();
	 // This shoudl be an instance of MoleculizerRulesFile

	 pArgs = PyTuple_New(1);
	 PyTuple_SetItem(pArgs, 0, pFileName);

	 pValue = PyObject_CallObject(pFunc, pArgs);
	 PyErr_Print();
	 
	 PyObject_CallMethod(pValue, "DEBUGPRINT", NULL);
	 PyErr_Print();

	 
       }
     else
       {
	 std::cout << "Nope!!!" << std::endl;
       }


     Py_DECREF( pName );
     Py_DECREF( pFileName );
     Py_DECREF( pArgs );
     Py_DECREF( pValue );


//       {
// 	pValue = PyObject_CallObject( pFunc, pValue);
// 	std::cout << "Tan of 100 is " << PyFloat_As( pValue ) << std::endl;

// 	// Py_DECREF(pValue);
	
//       }
//     else
//       {
// 	std::cout << "Error!" << std::endl;
// 	PyErr_Print();
//       }
    

    // Py_DECREF( pModule );
    // Py_DECREF( pName);

    //    Py_Finalize();
    
    return;
    
  }


  PythonRulesManager::PythonRulesManager()
    :
    _isInitialized( true )
  {

    std::cout << "Hello from the Python Rules Manager" << std::endl;
    std::cout << "Begin executing python experiments in the constructor()..." << std::endl;
    std::cout << "##############################" << std::endl;



    this->DEBUG_doInterestingStuff();


    std::cout << "##############################" << std::endl;

  }

  std::string
  PythonRulesManager::getXmlString() const
  {
    return std::string("");
  }



  void PythonRulesManager::addParameterStatement(const std::string& paramLine)
  {
    if (!isInitialized() ) throw std::exception();

    char functionName[] = "addParameterStatement";
    PyObject_CallMethod(mzrFileConverterClassInst, functionName, "s", paramLine);
  }

  void PythonRulesManager::addModificationStatement(const std::string& modLine)
  {
    if (!isInitialized() ) throw std::exception();


    PyObject_CallMethod(mzrFileConverterClassInst, "addModificationStatement", "s", modLine);

  }

  void PythonRulesManager::addMolsStatement(const std::string& molsLIne)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addMolsStatement", "s", molsLIne);

  }

  void PythonRulesManager::addAllostericPlexStatement(const std::string& alloPlexLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addAllostericPlexStatement", "s", alloPlexLine);

  }

  void PythonRulesManager::addAllostericOmniStatement(const std::string& alloOmniLine)
  {

    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addAllostericOmniStatement", "s", alloOmniLine);

  }

  void PythonRulesManager::addDimerizationGenStatement(const std::string& dimerGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addDimerGenStatement", "s", dimerGenLine);
  }

  void PythonRulesManager::addOmniGenStatement(const std::string& omniGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addOmniGenStatement", "s", omniGenLine);
  }

  void PythonRulesManager::addUniMolGenStatement(const std::string& uniMolGenLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addUniMolGenLine", "s", uniMolGenLine);
  }

  void PythonRulesManager::addSpeciesStreamStatement(const std::string& speciesStreamLine)
  {
    if (!isInitialized() ) throw std::exception();

    PyObject_CallMethod(mzrFileConverterClassInst, "addSpeciesStreamLine", "s", speciesStreamLine );
  }


}
