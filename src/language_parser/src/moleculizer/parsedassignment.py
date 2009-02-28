from parsedmzrtoken import ParsedMzrToken

class ParsedAssignment( ParsedMzrToken ):
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

    @staticmethod
    def AssignmentPassesSanityCheck( token):
        return "=" in token and not "->" in token and not "<-" in token

    def __init__(self, assignmentLine):
        ParsedMzrToken.__init__(self, assignmentLine)
        
        self.lhs = 0


        self.has_value = False
        self.has_expression = False

        self.rhs_value = 0
        self.rhs_expression = 0

        self.parse()

    def parse(self):
        token = self.getOriginalLine()
        
        if not self.AssignmentPassesSanityCheck(token):
            raise Exception()

        self.lhs = token.split("=")[0]

        rhsToken = token.split("=")[1]

        try:
            number = float( rhsToken )
            self.has_value = True
            self.rhs_value = number
        except:
            self.has_expression = True
            self.rhs_expression = rhsToken


    def getName(self):
        return self.lhs

    def getValue(self):
        # For now, don't worry about expressions.
        if not self.has_value:
            raise Exception()
        
        return self.rhs_value

    def isAssignment(self):
        return True
