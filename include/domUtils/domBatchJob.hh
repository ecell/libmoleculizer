/////////////////////////////////////////////////////////////////////////////
// Moleculizer - a stochastic simulator for cellular chemistry.
// Copyright (C) 2001  Walter Lawrence (Larry) Lok.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//    
// Contact information:
//   Larry Lok, Research Fellow          Voice: 510-981-8740
//   The Molecular Sciences Institute      Fax: 510-647-0699
//   2168 Shattuck Ave.                  Email: lok@molsci.org
//   Berkeley, CA 94704
/////////////////////////////////////////////////////////////////////////////

#ifndef DOMBATCHJOB_H
#define DOMBATCHJOB_H

#include <iostream>
#include "domUtils/domUtils.hh"

namespace domUtils
{
  class domBatchApp
  {
  public:
    virtual
    ~domBatchApp(void)
    {}

    virtual int
    run(void) throw(std::exception) = 0;
  };
    
  template<class batchAppClass>
  class domBatchJob
  {
    int argCount;
    char** argVector;
    
  public:
    domBatchJob(int argc, char** argv) :
      argCount(argc),
      argVector(argv)
    {}

    int
    parseNrun(void)
    {
      try
	{
	  // Create a non-validating parser.
	  xmlpp::DomParser parser;
	  parser.set_validate(false);

	  parser.parse_stream(std::cin);

	  // Xmlpp doc says this test is to "see if a document has been
	  // parsed."
	  if(parser)
	    {
	      // Create the application from the input document.
	      batchAppClass theApp(argCount,
				   argVector,
				   parser.get_document());

	      // Run the application.
	      return theApp.run();
	    }
	  else
	    throw noDocumentParsedXcpt();
	}
      catch(const std::exception& rExcept)
	{
	  std::cerr << rExcept.what()
		    << std::endl;
	  return 1;
	}
      return 0;
    }
  };

  template<class batchAppClass>
  class twoDomJob
  {
    int argCount;
    char** argVector;

  public:
    twoDomJob(int argc, char** argv) :
      argCount(argc),
      argVector(argv)
    {}

    int
    parseNrun(void)
    {
      try
	{
	  // Create a non-validating parser for standard input.
	  // I'm not entirely sure that one can't reuse a parser.
	  // Or maybe should.
	  xmlpp::DomParser stdInParser;
	  // Looking at the headers, I probably don't need to do
	  // this.
	  stdInParser.set_validate(false);

	  // Parse standard input.
	  stdInParser.parse_stream(std::cin);

	  // Create a non-validating parser for the auxiliary file.
	  xmlpp::DomParser auxFileParser;
	  auxFileParser.set_validate(false);

	  // Get the name of the auxiliary file from the command line.
	  if(argCount < 2) throw(insufficientArgsXcpt(argCount,
						      2));
	  const std::string auxFileName(argVector[1]);
	      
	  // Parse the file.
	  auxFileParser.parse_file(auxFileName);

	  if(stdInParser && auxFileParser)
	    {
	      // Create the application from the two input documents.
	      batchAppClass theApp(argCount - 1,
				   argVector + 1,
				   stdInParser.get_document(),
				   auxFileParser.get_document());

	      // Run the application.
	      return theApp.run();
	    }
	  else
	    throw noDocumentParsedXcpt();
	}
      catch(const std::exception& rExcept)
	{
	  std::cerr << rExcept.what()
		    << std::endl;
	  return 1;
	}
      // Technically I don't think I need this, but I wonder if gcc
      // could tell....
      return 0;
    }
  };
}

#endif // DOMBATCHJOB_H
