#ifndef __PYTHON_RULES_MANAGER_HH
#define __PYTHON_RULES_MANAGER_HH

#include <Python.h>
#include <string>


namespace mzr
{

  class PythonRulesManager
  {
  public:
    PythonRulesManager();
    ~PythonRulesManager();

    std::string getXmlString() const;

    std::string 
    addRulesFile( const std::string& rulesFile);
    
    std::string 
    addRulesString( const std::string& rulesString);

    void addParameterStatement(const std::string& paramLine);
    void addModificationStatement(const std::string& modLine);
    void addMolsStatement(const std::string& molsLIne);
    void addAllostericPlexStatement(const std::string& alloPlexLine);
    void addAllostericOmniStatement(const std::string& alloOmniLine);
    void addDimerizationGenStatement(const std::string& dimerGenLine);
    void addOmniGenStatement(const std::string& omniGenLine);
    void addUniMolGenStatement(const std::string& uniMolGenLine);
    void addSpeciesStreamStatement(const std::string& speciesStreamLine);

    bool isInitialized() const
    {
      return _isInitialized;
    }

  private:

    bool _isInitialized;

    char* getFileStringFunctionName;
    char* addWholeRulesFileFunctionName;
    char* addWholeRulesStringFunctionName;

    char* paramPythonFunctionName;
    char* modPythonFunctionName;
    char* molPythonFunctionName;
    char* alloPlexPythonFunctionName;
    char* alloOmniPythonFunctionName;
    char* dimerGenPythonFunctionName;
    char* omniGenPythonFunctionName;
    char* uniMolGenPythonFunctionName;
    char* speciesStreamPythonFunctionName;
    char* explicitSpeciesPythonFunctionName;
    char* singleStringTuple;

    PyObject* mzrFileConverterClassInst;

    void DEBUG_doInterestingStuff();

  };


}


#endif
