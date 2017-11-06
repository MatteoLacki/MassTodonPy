from operator import itemgetter
from collections import namedtuple

test_dict = {"a": 100, "b": 12, "c": 13}

for key, val in sorted(test_dict.items(), key=itemgetter(1)):
    print(key, val)


Point = namedtuple('Point', ['x', 'y'])
p = Point(x=10, y="b")
p.x
p.y

p
x, y = p
