def mytree(r):
	return [r,[],[]]
def insterLef(r,L):
	return [r,[mytree(L)],[]]
def insterRig(r,R):
	return [r,[],[mytree(R)]]











#file = open("a_sample_tree")
#for line in file:
#    pass print(line)
