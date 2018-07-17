class T(object):
    def __init__(self):
        self.diff = 1.0

    def iterate(self):
        i = 0
        while True:
            yield i
            i += self.diff

    def test(self):
        for i in self.iterate():
            self.diff = 2*self.diff
            yield i

t = T()

i = t.test()
next(i)