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

#ifndef MZRXCPT_H
#define MZRXCPT_H

#include <iostream>
#include <sstream>
#include <exception>

namespace mzr
{
  // Essentially a prefix for messages identifying all the remaining
  // diagnostic traps, where there is almost nothing to say that would be
  // useful to the user.
  std::string
  internalXcptMsg(void);

  class mzrXcpt : public std::exception
  {
    std::string message;

  public:

    mzrXcpt(const std::string& rMessage) throw() :
      message(rMessage)
    {}

    mzrXcpt(const char* pMessage) throw() :
      message(pMessage)
    {}

    virtual ~mzrXcpt(void) throw()
    {}

    const std::string&
    getMessage(void) const throw()
    {
      return message;
    }

    // Reimplementation of std::exception::what.
    const char*
    what(void) const throw()
    {
      return message.c_str();
    }
  };

  class mzrWarning
  {
    std::string message;
  public:

    mzrWarning(const std::string& rMessage) throw() :
      message(rMessage)
    {}

    mzrWarning(const char* pMessage) throw() :
      message(pMessage)
    {}

    virtual ~mzrWarning(void) throw()
    {}

    const std::string&
    getMessage(void) const throw()
    {
      return message;
    }

    void
    issue(void) const throw()
    {
      std::cerr << "Warning: "
		<< getMessage()
		<< std::endl;
    }
  };

  class mzrMessage
  {
    std::string message;
  public:

    mzrMessage(const std::string& rMessage) throw() :
      message(rMessage)
    {}

    mzrMessage(const char* pMessage) throw() :
      message(pMessage)
    {}

    virtual ~mzrMessage(void) throw()
    {}

    const std::string&
    getMessage(void) const throw()
    {
      return message;
    }

    void
    issue(void) const throw()
    {
      std::cerr << getMessage()
		<< std::endl;
    }
  };

  class badDumpFileXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(const std::string& rBadFileName)
    {
      std::ostringstream msgStream;
      msgStream << "Could not open dump file `"
		<< rBadFileName
		<< "' for writing.";
      return msgStream.str();
    }
  public:
    badDumpFileXcpt(const std::string& rBadFileName) :
      mzrXcpt(mkMsg(rBadFileName))
    {}
  };

  class timeLimitExpiredXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(double clockTimeLimit,
	  double curSimTime)
    {
      std::ostringstream msgStream;
      msgStream << "Clock time limit of "
		<< clockTimeLimit
		<< " seconds exceeded at simulation time "
		<< curSimTime
		<< ".";
      return msgStream.str();
    }
  public:
    timeLimitExpiredXcpt(double clockTimeLimit,
			 double curSimTime) :
      mzrXcpt(mkMsg(clockTimeLimit,
		    curSimTime))
    {}
  };

  class eventQueueEmptyXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(double curSimTime)
    {
      std::ostringstream msgStream;
      msgStream << "Event queue empty at simulation time "
		<< curSimTime
		<< ".";
      return msgStream.str();
    }
  public:
    eventQueueEmptyXcpt(double curSimTime) :
      mzrXcpt(mkMsg(curSimTime))
    {}
  };

  class nowIsNeverXcpt : public mzrXcpt
  {
  public:
    nowIsNeverXcpt(void) :
      mzrXcpt("Earliest event is at 'infinite' simulation time.")
    {}
  };

  class negativePopXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(const std::string& rSpeciesName)
    {
      std::ostringstream msgStream;
      msgStream << "Population of species "
		<< rSpeciesName
		<< " became negative.";
      return msgStream.str();
    }
  public:
    negativePopXcpt(const std::string& rSpeciesName) :
      mzrXcpt(mkMsg(rSpeciesName))
    {}
  };

  // Thrown in featureMap.hh, featureMap<...>::addFeature.
  class featureAlreadyMappedXcpt : public mzrXcpt
  {
    static std::string
    mkMsg(void)
    {
      std::ostringstream msgStream;
      msgStream << mzr::internalXcptMsg()
		<< "featureSpec is already mapped.";
      return msgStream.str();
    }
  public:
    featureAlreadyMappedXcpt(void) :
      mzr::mzrXcpt(mkMsg())
    {}
  };
}

#endif // MZRXCPT_H
