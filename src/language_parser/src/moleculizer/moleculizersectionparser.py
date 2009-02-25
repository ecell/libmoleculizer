###############################################################################
# Libmoleculizer
# Copyright (C) 2007, 2008, 2009 The Molecular Sciences Institute
#
# Moleculizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Moleculizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Moleculizer; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Original Author:
#   Nathan Addy, Scientific Programmer 
#   The Molecular Sciences Institute    Email: nathan.addy@gmail.com
#   
###############################################################################

class MoleculizerSection( object ) :
    HAVE_EE = False

    @classmethod
    def getExpressionEvaluator(cls):
        return cls.HAVE_EE

    @classmethod
    def setExpressionEvaluator( cls, expression_evaluator):
        cls.HAVE_EE = expression_evaluator
    
    @classmethod
    def getValue( cls, expression ):
        if (HAVE_EE):
            value = self.evaluateExpression( expression )
            return value
        else:
            MissingExpressionEvaluatorException()
                
    def parse(self):
        raise UnimplementedParceFunction( "Unimplemented parse function in class '%s'" % str( self.__class__ ).split("'")[1] )
