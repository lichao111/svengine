#!user/bin/env python

import tempfile,os,string

generation=3 #add a warning not let genaration bigger than real generation

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
#def trefil2varfile()

def P(Node):	#find Node's parent
	i=Node.i-1
	j=Node.j/2
	return globals()['N'+str(i)+str(j)]
def anchor(Node,n):	#find Node's anchor
	for i in range(n):
		Node=P(Node)
	return Node

def C(Node):	#find Node's children
	i=Node.i+1
	j1=Node.j+1
	j2=Node.j+2	
	return [locals()['N'+str(i)+str(j1)],locals()['N'+str(i)+str(j2)]]

def depth(j,generation):	#the depth of j brunch, generation minu the number of None in this brunch
	assert j in range(2**generation),	'%d  brunch in all' %(2**generation)	
#	print Nonelist	
	i=generation
	depth=generation
#	for l in Nonelist:
#		if  not in l[1]:
#			return generation/2
#		else:
			
	while [i,j] in Nonelist:
#		print [i,j]
		i=i-1
		j=int(j/2)
		depth=depth-1
	return depth		

def get_idname (d,Node):	#d=locals(),but not here
	for key in d:
		if d[key]==Node:
			return key
		
	
Nodelist=[]	#store all Node class
id_Nodelist=[]	#store all Node id, for example N00,N10...
Nonelist=[]	#store the 'None' Node
VAR=[]	#store all kinds of var name
with open('/home/lichao/Documents/github/svengine/exam_tree') as treefile:
#	treefile=file

#	with open('/home/lichao/Documents/github/svengine/tree1','w+') as tree1:
#tree1=open('/home/lichao/Documents/github/svengine/tree1','w+')
#tree1=tempfile.TemporaryFile()
#tree1=open('tree1','w+')
		for g,line in enumerate(treefile):
			if '/' not in line and g<=2*generation:
				#tree1.write(line)
				line=line.expandtabs().strip().split(' ')
				while '' in line:
					line.remove('')		
				l=len(line)
				
				if g==0:# the root node
					print line	
					assert  line==['N'],	"the tree should use N as root Node"
					locals()['N'+'0'+'0']=Node(0,0,{'Normal':None},1)
					Nodelist.append(Node(0,0,{'Normal':None},1))
				else:
					for j in range(l) :		
						
						if line[j]=='None':	#check the 'None' condition
							#print line[j],type(line[j])							
							Nonelist.append([g/2,j])
#							print Nonelist							
							locals()['N'+str(g/2)+str(j)]=Node(g/2,j,locals()['N'+str(g/2-1)+str(int(j/2))].v,0.5)
							Nodelist.append(Node(g/2,j,locals()['N'+str(g/2-1)+str(int(j/2))].v,0.5))
						else:
						
							
						
							

							line[j]=line[j].split(';')
							
							var=line[j][0].split(',')
							freq=line[j][1]
							#print	VAR.append			
						
							#print var.split(':')[0],var.split(':')[1],var,type(var)
						
							if var[0][0]=='N' and var[0][1].isdigit() :
								assert len(var)==1,	"Normal Node only have one type 'N'"	
								locals()['N'+str(g/2)+str(j)]=Node(g/2,j,{'Normal':None},freq)	
								Nodelist.append(Node(g/2,j,{'Normal':None},freq))
							else:
								for i in range(len(var)):
									if var[i] not in VAR:
										VAR.append(var[i])#f  var not in VAR:
						
	#						line[j]=line[j].split(';')
	#						b=line[j][0]
								varlist={}
							#freq=line[j][1]
							#print varlist
								
								map(lambda x:varlist.setdefault(x.split(':')[0], int(x.split(':')[1])), var)
						
								locals()['N'+str(g/2)+str(j)]=Node(g/2,j,varlist,freq)
								Nodelist.append(Node(g/2,j,varlist,freq))
					#if line[j].startwith('N'):
					#	locals()['N'+g/2+j]=Node(g/2,j,varlist,freq)
#print Nodelist,	'this the all Node!'
assert len(Nodelist)==2**(generation+1)-1,	'the number of Node should equal to 2**(generation+1)-1, the treefile may include error.'	
print VAR
#for i in range(2**generation):
#	print locals()['N'+str(generation)+str(i)].v
#for i in range(15):
#	print Nodelist[i].i,Nodelist[i].j,Nodelist[i].v,Nodelist[i].a	
#print N32.v,N21.v,N32.a
#print locals()		
#print get_name(Node(0,0,{'Normal':None},1)),get_name(N00)#

#print Node(0,0,{'Normal':None},1),N00	#why False
#print locals()['N00'] is Node(0,0,{'Normal':None},1)
#print id(N00),'\t',id( Node(0,0,{'Normal':None},1))
#print N00==Node(0,0,{'Normal':None},1)
for i in range(generation+1):
	for j in range(2**i):
		id_Nodelist.append(get_idname(locals(),locals()['N'+str(i)+str(j)]))
print id_Nodelist
#print locals(),'\n',globals()
#print locals()['N'+str(g/2)+'1']
#print get_idname(locals(),anchor(P(locals()['N'+str(g/2)+'1'],3)
#print anchor(locals()['N'+str(g/2)+'1'],3)
for i in range(2**generation):
	print depth(i,generation)

def freq(vari,generation):		#calculte the allel freq
	cell_all=0	#the number of all cell
	cell_vari=0	#the number of cell inculde vari: if the vari is heterozygote,the number be half
	count_brunch=[]	#the number of cell each brunch
	brunch_vari=[]	#the brunch id which include vari 
	for i in range(2**generation):
		a=1
		for j in range (generation):
			
			a=a*(anchor(globals()['N'+str(generation)+str(i)],j).a)
#			print a
		
		count_brunch.append(2**(depth(i,generation))*a)
#	print count_brunch	
	cell_all =round(sum(count_brunch),2)
	haplist=[]
	for i in range(2**generation):
		if vari in globals()['N'+str(generation)+str(i)].v.keys():
			haplist.append(globals()['N'+str(generation)+str(i)].v[vari])
			brunch_vari.append(i)
	if haplist.count(haplist[0])!=len(haplist):	
		print "Warning! The haplotype of each var must keep unchange,you better check your %s  "	%vari
	hap=haplist[0]
	assert hap in [0,1],	"the haplotype must 1(homozygote) or 0(heterozygote)"
	for j in brunch_vari:
		cell_vari=(cell_vari+count_brunch[j])	#0-->0.5,1-->1
	cell_vari=round(cell_vari*(float(hap+1)/2),2)
								
	return cell_vari,cell_all,round(cell_vari/cell_all,2)				
		
#print freq('var1')#,freq('var2'),freq('var3'),freq('var4'),freq('var5')
def VAR_frq(VAR,generation):
	VAR_frq={}
	for var in VAR:
		VAR_frq[var.split(':')[0]]=freq(var.split(':')[0],generation)	
	return VAR_frq
print VAR_frq(VAR,generation)	
	

#print id_Nodelist[1].a

#Nodelist1=[]
#print N11.i,N11.j,N11.v,N11.a				
#for g in range(generation+1):
#	for j in range(2**g):
#		Nodelist1.append(locals()['N'+str(g)+str(j)])
#print [N00,N10,N11]
#print ''
#for line in tree1:
#				print g/2,line  #g/2 because '/' was remove
#quit()
#		with open('/home/lichao/Documents/github/svengine/tree2','w+') as tree2:
#			for g,line in enumerate(tree1):
#				print g,line
#				if g<=generation :
#					tree2.write(line)
#
#	quit()
	#tree2=open('/home/lichao/Documents/github/svengine/tree2')


#A00=Node(0,0,{'var1':1},0.5)
#print A00,A00.i,A00.j,A00.v,A00.a,type(A00.v)	
			#for g,line in enumerate(tree2):
	#line=line.expandtabs().strip().split(' ')
	#string.translate(line, del=" ")
	#line=line.strip().split(' ')
	#line=line.split(' ')
	
#	while '' in line:
#		line.remove('')
#	l=len(line)
#	for j in range(l):
#		line[j]=line[j].split(';')
#		varlist=line[j][0]
#		print varlist
#		a={}
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
