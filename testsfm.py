

import statistics as stats

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg


class class_B(object):
    def __init__(self,name="Alex"):
        self.name = name
    def say_my_name(self):
        print self.name

class class_A(class_B):
    def __init__(self,name,code=100):
        class_B.__init__(self,name)
        self.code = code
    def say_goodbye(self):
        print "goodbye"

if __name__ == "__main__":

    A = class_A(name="Sean")
    A.say_my_name()
    A.say_goodbye()
    print A.code