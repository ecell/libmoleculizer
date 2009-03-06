#include <Python.h>

#include "pythonRulesManager.hh"


#include <iostream>

namespace mzr
{

  void PythonRulesManager::DEBUG_doInterestingStuff()
  {
    PyObject *pName, *pModule, *pDict, *pFunc, *pValue;

    Py_Initialize();

    pName = PyString_FromString( "math" );
    pModule = PyImport_Import( pName );
    pDict = PyModule_GetDict( pModule );
    pFunc = PyDict_GetItemString( pDict, "pi" );
    std::cout << "Tan of 100 is " << PyFloat_AsDouble( pFunc ) << std::endl;

//     if (PyCallable_Check( pFunc ) )
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
    
	PyErr_Print();
    // Py_DECREF( pModule );
    // Py_DECREF( pName);

    //    Py_Finalize();
    
    return;
    
  }


  PythonRulesManager::PythonRulesManager()
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
  {}

  void PythonRulesManager::addModificationStatement(const std::string& modLine)
  {}

  void PythonRulesManager::addMolsStatement(const std::string& molsLIne)
  {}

  void PythonRulesManager::addAllostericPlexStatement(const std::string& alloPlexLine)
  {}

  void PythonRulesManager::addAllostericOmniStatement(const std::string& alloOmniLine)
  {}

  void PythonRulesManager::addDimerizationGenStatement(const std::string& dimerGenLine)
  {}

  void PythonRulesManager::addOmniGenStatement(const std::string& omniGenLine)
  {}

  void PythonRulesManager::addUniMolGenStatement(const std::string& uniMolGenLine)
  {}

  void PythonRulesManager::addSpeciesStreamStatement(const std::string& speciesStreamLine)
  {}


}
