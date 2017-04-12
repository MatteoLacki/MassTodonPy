

class TT:
    def __init__(self):
        self.arg = 10

    def m(self, arg=None):
        if arg is None:
            arg = self.arg
        return arg**2

tt = TT()

tt.m()
tt.m(3)
