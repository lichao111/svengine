#!user/bin/env python
#Representing a Tree Using a List, Binary Search Tree and Inorder Traversal
import sys

def treeInsert(s,a):
    if s == []: # first time through
        s.append(a)
        s.append([])
        s.append([])
    else:
        while True:
            if a == s[0]: # found
                return 1
            elif a < s[0]: # left
                if s[1] == []:
                    # not found, insert
                    s[1] = [a,[],[]]
                    return 0
                else:
                    s = s[1] # go left
            elif a > s[0]: # right
                if s[2] == []:
                    # not found, insert
                    s[2] = [a,[],[]]
                    return 0
                else:
                    s = s[2] # go right

def binarySearchTree(m):
    r = [] # build tree from scratch
    for i in range(0,len(m)):
        stat = treeInsert(r,m[i])
    return r

def inorderTraversal(s):
    if (s[1] != []): # left subtree
        inorderTraversal(s[1])
    sys.stdout.write(str(s[0]) + " ")
    if (s[2] != []): # right subtree
        inorderTraversal(s[2])

def process(m):
    sys.stdout.write("INPUT DATA:        " + str(m) + '\n')
    r = binarySearchTree(m)
    sys.stdout.write("TREE AS STRING:    " + str(r) + '\n')
    sys.stdout.write("INORDER TRAVERSAL: ")
    inorderTraversal(r)
    sys.stdout.write("\n\n")

m = ["Mon","Tue","Wed","Thu","Fri","Sat","Sun",
     "Jan","Feb","Mar","Apr","May","Jun",
     "Jul","Aug","Sep","Oct","Nov","Dec",]
process(m)

n = [10.5,  3, 11.25, 17, 6.75, 7,  19,
     23.5, 14, 19.75, 12, 4,    1, 5.5]
process(n)
b=[10.5,3,11,11,17]
process(b)

