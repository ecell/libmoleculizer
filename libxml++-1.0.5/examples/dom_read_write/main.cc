// -*- C++ -*-

/* main.cc
 *
 * Copyright (C) 2002 The libxml++ development team
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <libxml++/libxml++.h>

#include <iostream>


int
main(int argc, char* argv[])
{
  //Parse command-line arguments:
  std::string filepath_in, filepath_out;
  if(argc > 1 )
    filepath_in = argv[1]; //Allow the user to specify a different XML file to parse.
  else
    filepath_in = "example.xml";

  if(argc > 2 )
    filepath_out = argv[2]; //Allow the user to specify a different output file.
  else
    filepath_out = "example_output.xml";
        
  try
  {
    xmlpp::DomParser parser;
    parser.set_validate();
    parser.parse_file(filepath_in);
    if(parser)
    {
      //Write it out again.
      xmlpp::Document* document = parser.get_document();
      if(document)
        document->write_to_file(filepath_out);
    }
  }
  catch(const std::exception& ex)
  {
    std::cout << "Exception caught: " << ex.what() << std::endl;
  }

  return 0;
}

