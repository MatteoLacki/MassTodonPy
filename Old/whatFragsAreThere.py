 len(fasta)
 frags = set()
 for re,e,s in res:
     for r in re:
         frags.add(r['molType'])

 'z76' in frags
 Forms = makeFormulas(fasta=fasta, Q=Q, fragType='cz')
 list(Forms.makeMolecules())
 frags_made = set()
 for t,_,b,_,_ in Forms.makeMolecules():
     frags_made.add(t[0]+str(b))
 'z76' in frags_made
 res
