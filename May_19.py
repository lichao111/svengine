#!user/bin/env python

import tempfile,os,string

generation=3  #add a warning not let genaration bigger than real generation

#tree1=open('/home/lichao/Documents/github/svengine/tree1','w+b')
#for G,line in enumerate(o):
#	print G,line
#home=tempfile.TempFile()
#print 

class Node:
	
	def __init__(self,i,j,v,a):#i:int,j:int,v:dict,a:range[0,1]
		self.i=int(i)
		self.j=int(j)
		self.v=v
		self.a=float(a)
#treefile=file('/home/lichao/Documents/github/svengine/exam_tree')
with open('/home/lichao/Documents/github/svengine/exam_tree') as file:
	treefile=file

	with open('/home/lichao/Documents/github/svengine/tree1','w+') as tree1:
#tree1=open('/home/lichao/Documents/github/svengine/tree1','w+')
#tree1=tempfile.TemporaryFile()
#tree1=open('tree1','w+')
		for g,line in enumerate(treefile):
			if '/' not in line and g<=2*generation:
				tree1.write(line)
				line=line.expandtabs().strip().split(' ')
				while '' in line:
					line.remove('')		
				l=len(line)
				for j in range(1,l) :#check 'N'	 condition			
					line[j]=line[j].split(';')
					varlist=line[j][0]
					a={}
					#freq=line[j][1]
					#print varlist
					map(lambda x:a.setdefault(x.split(':')[0], x.split(':')[1]), varlist)
					
					quit()
					#if line[j].startwith('N'):
					#	locals()['N'+g/2+j]=Node(g/2,j,varlist,freq)
		
	
				
					


#for line in tree1:
#				print g/2,line  #g/2 because '/' was remove
#quit()
#		with open('/home/lichao/Documents/github/svengine/tree2','w+') as tree2:
#			for g,line in enumerate(tree1):
#				print g,line
#				if g<=generation :
#					tree2.write(line)
#
	quit()
	#tree2=open('/home/lichao/Documents/github/svengine/tree2')


#A00=Node(0,0,{'var1':1},0.5)
#print A00,A00.i,A00.j,A00.v,A00.a,type(A00.v)	
			#for g,line in enumerate(tree2):
	#line=line.expandtabs().strip().split(' ')
	#string.translate(line, del=" ")
	#line=line.strip().split(' ')
	#line=line.split(' ')
	
	while '' in line:
		line.remove('')
	l=len(line)
	for j in range(l):
		line[j]=line[j].split(';')
		varlist=line[j][0]
		print varlist
		a={}
	#	map(lambda x:a.setdefault(x.split(':')[0], x.split(':')[1]), varlist)
	#	if line[j].startwith('N'):
	#		locals()['N'+g+j]=	
		
	
		
	#print line	
		#if :
			
	#print line,l
#for g,line in enumerate(tree2):
#	print g,line
#for line in tree1:
#	print line	
	
#tree1.close()
