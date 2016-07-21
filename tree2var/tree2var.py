#!/usr/bin/env python

import random,copy,csv,collections,time,argparse,os

start_time=time.time()

varfile_names = ['#VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)']

def find_event(t,event=[]):
	if t:
		if t[0] not in event:
			event.append(t[0])
		find_event(t[1],event)
		find_event(t[2],event)
#	else:
	return event
	
#def get_fre(t,s=0):
#	global m,n
#	if t[1] or t[2] :
#		if s==0:
#			fre[0]=random.uniform(0,1)
#			fre[1]=random.uniform(0,1)
#			fre[2]=random.uniform(0,1)
#		if s==1:
#			fre[m]=random.uniform(0,1)
#			fre[m*10+1]=random.uniform(0,1)
#			fre[m*10+2]=random.uniform(0,1)
#			m=m*10+1
#		if s==2:
#			fre[n]=random.uniform(0,1)
#			fre[n*10+1]=random.uniform(0,1)
#			fre[n*10+2]=random.uniform(0,1)		
#			n=n*10+1
#		
#		get_fre(t[1],1)
#		get_fre(t[2],2)
#	return fre

def is_emp(tree,t):
	
	w=len(str(t))
	for i in range(w):
		a=int(str(t)[::-1][w-1])
#		print a
		tree = tree[a]
#		print tree
		w = w-1
	return tree[1] == []	
	
def g(a):	#get the children node number
	b=copy.deepcopy(a)
	for i in range(len(a)):
		b.extend([b[i]*10+1,b[i]*10+2])
	b=list(set(b))
#	print b
#	print a
	for i in a:
		b.remove(i)
	return b
	
def tai_tre(tree):	#find the leaf node ,default the first floor [1,2] not empty
	num_index=0;
	t=[1,2]
	b=[]
	tail=[]
	T=[1,2]
	while 1:
		for i in t:	
			try:
				if is_emp(tree,i):
					tail.append(i)
				#tail=list(set(tail))
#				print 'tail=',tail
			except IndexError:
#				print t,'i',i
				T.remove(i)
#				print "T=",T
				#num_index += 1
		#	num_index += 1
		#if num_index == len(t):
		#if num_index == len(t):
#		print "t=",t
		if T==[] :
			return tail#break
 		t=g(T)
 		T=copy.copy(t)

def mate(tree,event,tail):	#mate the node with the event ##find the min node 
	node=[0]
	for t in tail:
		for i in range(len(str(t))):
			while t:
				node.append(t)
				t = t/10
	node=list(set(node))			
#	print 'node',node	
	mat={}
	for i in event:
#		print i,type(i) 	
		for j in node:
#			print '--------',j,'---------'
			tree1=copy.deepcopy(tree)
#			print tree1
			w=len(str(j))
			for k in range(w):
				a=int(str(j)[::-1][w-1])
#				print a
				tree1 = tree1[a]
#				print tree
				w = w-1
#			return tree[0] == i
#			print tree1,type(tree1[0]),tree1[0],i
			if tree1[0] == i:
#				print 'YES'
				mat[j]=i
	T={}	
	for i in event:
		T[i]=[]		
	for i in mat:
		for j in event:
			#T[j]=[]
			#print mate[i],
			if mat[i] == j:
				T[j].append(i)
	for i in T:
		T[i]=min(T[i])
#	print 'T',T	
	T=dict((value,key) for key,value in T.iteritems())
	node=[i for i in T.keys()]			
	return [T,node]

def che_odd(a): #check a number is odd or even
	if a%2 ==0:
		return False
	else:
		return True

def get_fre(mat,tail):	#calculte the frequence of the event
	fre={}
	fre[0]=round(random.uniform(0,1),2)
	for i in tail:
 
#		print fre
		w = len(str(i))
		for  j in range(w-1):
			#print i/(10**(w-1))
			fre[i/(10**(w-1))]=round(random.uniform(0,1),2)
			w -= 1
			#fre[i] *= fre[i/(10**(w-1))]
#	print 'fre',fre
	fre2={}
	for i in tail:
		fre2[i]=1
		l=copy.copy(i)
		
#		if che_odd(i):
#			fre[i] = 1-fre[0]
#		else:
#			fre[i] = fre[0]
#		print l
		w=len(str(i))
		for k in range(w):
			#l=l/(10**(w-1))
			if i == 0:
				fre2[i]=fre[0]
				
			elif che_odd(l):
				fre2[i] *= 1-fre[l/(10)]
				fre2[i] = round( fre2[i],2)
#				print 'fre2[%s]' %(i),'*=','1-fre[%s]'%(l/10)
			else:
				fre2[i] *= fre[l/(10)]
				fre2[i] = round( fre2[i],2)
#				print 'fre2[%s]'%i,'*=','fre[%s]'%(l/10)
			l = l/10
#	print 'fre2',fre2		
	fre1={}
	for i in mat.keys():
		for j in fre2.keys():
			if j == i:
				fre1[mat[i]]=fre2[j]
			else:
				pass 
			
#			fre1[i]=	
#	print "fre",fre
	return fre1

def var_fre(varfile,varfre,commentstring='#'):
	varlist=collections.OrderedDict()
	for line in varfile:
		if not line.startswith(commentstring):
			if line.endswith('\n'):
				line=line[::-1][1:][::-1]
			line = line.split('\t')
			varlist[line[0]]=line
	tes_var=[]		
	for i in varlist:
		for j in varfre:
			if i == j:
				tes_var.append(i)
				if varlist[i][2]== '0':
					varlist[i]=varlist[i][:2]+[varfre[j]]+varlist[i][2:]
				elif varlist[i][2]== '1':
					varfre[j]=str(float(varfre[j])/2)
					varlist[i]=varlist[i][:2]+[varfre[j]]+varlist[i][2:]

				else:
					print 'The Hap of the event may error,only in 1 or 0'
					quit()	
	assert tes_var == [ i for i in varlist], "please check if the event in tree math with the event in var file "
	return varlist

def main(args):
#	print 'OK'
	tree=eval(args.tree)
	varfile=args.varfile
	oprefix=args.oprefix
	tail=tai_tre(tree)

	event=find_event(tree)
	print 'total event:',event
	mat=mate(tree,event,tail)[0]
	node=mate(tree,event,tail)[1]
	fre1=get_fre(mat,node)
	varlist=var_fre(varfile,fre1)
	event_fre=copy.deepcopy(fre1)
	[event_fre.pop('N') if 'N' in event_fre else event_fre]
	print 'the frequency of each event:',event_fre
	#oname=varfile.name.rsplit(',')[0]
	if not oprefix:
		oprefix=varfile.name.rsplit(',')[0]
#		var_out=csv.write(open('%s'+'_out.var','w')) %(oname)
	var_out=csv.writer(open(oprefix+'_out.var','w'),delimiter='\t') 
	var_out.writerow(varfile_names)
#	print oprefix
	for var in varlist:
		var_out.writerow(varlist[var])

	#print type(tree),tree,type(varfile),type(oprefix)
	#print type(oprefix)

def run():
	parser = argparse.ArgumentParser(description="translate a binary tree to varfile can be accepted by fasforge and muteforge.")
	parser.add_argument('tree',metavar='tree', type=str,help="input the evolution tree, required")
	parser.add_argument('varfile', metavar='varfile', type=argparse.FileType('rU'),help="inputfile of variant specification,with no frequency,required")
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, help="prefix for output files")
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
	print "tree2var have done! Total runtime: ", time.time()-start_time, " seconds"
