nodes = [[1,2],[3],[5],[5],[6],[],[]]
cluster = [0]*len(nodes)


class Queue(object):
    def __init__(self, queue=None):
        if queue is None:
            self.queue = []
        else:
            self.queue = list(queue)
    def dequeue(self):
        return self.queue.pop(0)
    def enqueue(self, element):
        self.queue.append(element)


q = Queue()
q.enqueue(1)
q.enqueue(10)

print q.dequeue

