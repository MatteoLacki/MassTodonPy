from heapq import heappush, heappop

h = []

heappush(h, (1,'a'))
heappush(h, (0,'b'))
heappush(h, (10,'ba'))
heappush(h, (-10,'cba'))

print h

heappop(h)
print h