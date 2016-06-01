#user/bin/env python
import mutforge as mf
import collections
freq1 = [ float(f) for f in [4,5,7]]
print freq1
for f in freq1:
		assert f>=4 and f<=7, "mixfreq must be a float between 0 and 1"
print "OOOOOOOOOOOOOOOOOO"
args="4,5,6,7,8,"
b=args.split(",")
print args,type(args),type(b),type({4,5}),type([4,5])
print '++++++++++++++++++'
try:
	8==9
except no:
	raise Exception("No")
print "-__________________test______________"
varlist = collections.OrderedDict()
dic={'nligation':0,'varlist':varlist}#prefix':oprefix,'reffile':reffile,'nploid':nploid,'plansize':plansize,'varcnt':0,'edgein':edgein,'gaptab':gaptab,'varfile':varfile,'parfile':parfile,'parlist':parlist,'nprocs':nprocs,'burnin':burnin,'mergemax':mergemax}
varlist=mf.file2var(mf.Bunch(dic))

