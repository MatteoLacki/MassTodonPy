from collections import Counter

class LinearCounter(Counter):
    def __add__(self, other):
        '''Add counts from two LinearCounters.

        >>> LinearCounter('abbb') + LinearCounter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})

        >>> LinearCounter({'H':-4, 'Y':10}) + LinearCounter({'H=2','Z':3})
        LinearCounter({'H': -2, 'Y': 10, 'Z': 3})

        '''
        if not isinstance(other, LinearCounter):
            return NotImplemented
        result = LinearCounter()
        for elem, count in self.items():
            newcount = count + other[elem]
            result[elem] = newcount

        for elem, count in other.items():
            if elem not in self:
                result[elem] = count
        return result

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __iadd__(self, other):
        '''Add counts from two LinearCounters.
        '''
        if not isinstance(other, LinearCounter):
            return NotImplemented
        for elem in self.keys():
            self[elem] += other[elem]
        for elem, count in other.items():
            if elem not in self:
                result[elem] = count
        return self

    def __isub__(self, other):
        '''Add counts from two LinearCounters.
        '''
        if not isinstance(other, LinearCounter):
            return NotImplemented
        for elem in self.keys():
            self[elem] -= other[elem]
        for elem, count in other.items():
            if elem not in self:
                self[elem] = -count
        return self

    def __mul__(self, scalar):
        '''Multiplies the values stored in the LinearCounter by the scalar.'''
        result = LinearCounter()
        for elem, count in self.items():
            result[elem] = count*scalar
        return result

    def __rmul__(self, scalar):
        '''Multiplies the values stored in the LinearCounter by the scalar.'''
        result = LinearCounter()
        for elem, count in self.items():
            result[elem] = count*scalar
        return result

# x = LinearCounter({'H':10, 'O':5})
# y = LinearCounter({'H':10.2, 'O':5.1})
# z = LinearCounter({'H':10.3, 'O':5.112, 'M':24232.232})
#
#
# x += y
# x
# x -= z
# x
# sum([x,y,z])
