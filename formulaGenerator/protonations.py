def protonate(Q,frag):
    a, b, c = {
        'precursor' : (1,0,1),
        'c' : (0,-1,0),
        'z' : (0,0,1)
    }[frag]
    for q in range(1,Q+a):
        for g in range(b,Q-q+c):
            yield (q,g)
