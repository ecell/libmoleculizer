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

#ifndef UTL_DOMJOBIMPL_H
#define UTL_DOMJOBIMPL_H

#include <libxml++/libxml++.h>
#include "utl/noDocumentParsedXcpt.hh"
#include "utl/insuffArgsXcpt.hh"

namespace utl
{
  namespace dom
  {
    template<class batchAppClass>
    int
    domBatchJob<batchAppClass>::
    parseNrun(void)
    {
      try
	{
	  // Create a non-validating parser.
	  typename xmlpp::DomParser parser;
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

    template<class batchAppClass>
    int
    twoDomJob<batchAppClass>::
    parseNrun(void)
    {
      try
	{
	  // Create a non-validating parser for standard input.
	  // I'm not entirely sure that one can't reuse a parser.
	  // Or maybe should.
	  typename xmlpp::DomParser stdInParser;
	  // Looking at the headers, I probably don't need to do
	  // this.
	  stdInParser.set_validate(false);

	  // Parse standard input.
	  stdInParser.parse_stream(std::cin);

	  // Create a non-validating parser for the auxiliary file.
	  typename xmlpp::DomParser auxFileParser;
	  auxFileParser.set_validate(false);

	  // Get the name of the auxiliary file from the command line.
	  if(argCount < 2)
	    throw insuffArgsXcpt::counts(argCount,
					 2);

	  // Parse the file.
	  const typename std::string auxFileName(argVector[1]);
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
      catch(const typename std::exception& rExcept)
	{
	  std::cerr << rExcept.what()
		    << std::endl;
	  return 1;
	}
      // Technically, I think this point in the code is not reachable.
      return 0;
    }
  }
}

#endif // UTL_DOMJOBIMPL_H
