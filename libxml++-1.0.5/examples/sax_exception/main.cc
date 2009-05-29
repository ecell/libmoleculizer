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

/* This example is nearly the same than sax_parser, but demonstrate how exception
 * can be thrown from the callback method as long as they herit from xmlpp::exception
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "myparser.h"
#include <iostream>

int main(int argc, char* argv[])
{
  MySaxParser parser;
  
  try
  {
    parser.parse_file("example.xml");
  }
  catch(const std::exception& ex)
  {
    std::cout << "Exception caught: " << ex.what() << std::endl;
  }

  return 0;
}

