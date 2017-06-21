def change_key(seq, q, p, name):
    if name[0]=='c':
        name = 'c'+str(int(name[1:])-1)
    return name, q, p
