import numpy as np
import re
str1="dataset/krogan2006extended.txt"
str2="matrix/Network_krogan2006extended.txt"
str3="matrix/Attribute_krogan2006extended.txt"

file1=open(str1)
file2=open(str2,'w')
filetemp=open("krogan2006extended_node.txt",'w')
print "create topological network!"
node=[]
for j in file1:
    temp1=j.split('	')[0]
    temp2=j.split('	')[1].rstrip('\n')
    if temp1 not in node:
        node.append(temp1)
        filetemp.write(temp1)
        filetemp.write('\n')
    if temp2 not in node:
        node.append(temp2)
        filetemp.write(temp2)
        filetemp.write('\n')
file1.close()
filetemp.close()
print node

file1=open(str1)

l = [([0] * len(node)) for x in range(len(node))]
for i in file1:
    temp1=i.split('	')[0]
    temp2=i.split('	')[1].rstrip('\n')
    q=0
    p=0
    for n in node:
        if n==temp1:
            a=node.index(n)
            q=1
            break
    for m in node:
        if m==temp2:
            b=node.index(m)
            p=1
            break
    if(q==1)and(p==1):
        l[a][b]=1
        l[b][a]=1
for i in range(len(node)):
    for j in range(len(node)-1):
        file2.write(str(l[i][j]))
        file2.write(' ')
    file2.write(str(l[i][len(node)-1]))
    file2.write('\n')
file2.close()

print "create attributed network!"

go=[]
file=open("krogan2006extended_go_information.txt")
file4=open("krogan2006extended_go_information_temp.txt",'w')
for i in file:
    node_name=i.split(' ',1)[0]
    node_go=i.split(' ',1)[1].rstrip('\n').rstrip(' ')
    file4.write(node_name)
    file4.write(' ')
    node_go=re.split(" |:",node_go)
    for j in node_go:
        if(j!='GO'):
            file4.write(j)
            file4.write(' ')
    file4.write('\n')
    for j in node_go:
        if j not in go:
            go.append(j)
file.close()
file4.close()
go.remove('GO')
go.sort()

file=open("krogan2006extended_go_information_temp.txt")
gov=[]
for i in file:
    node_name = i.split(' ', 1)[0]
    node_go = i.split(' ', 1)[1].rstrip('\n').rstrip(' ')
    node_go=node_go.split(' ')
    one = {}
    one['node_name']=node_name
    one['node_go']=node_go
    gov.append(one)
file.close()

attr = [[0 for col in range(len(go))] for row in range(len(node))]
file3=open(str3,'w')
for i in node:
    for j in gov:
        if i==j['node_name']:
            a=node.index(i)
            for z in j['node_go']:
                for q in go:
                    if z==q:
                        b=go.index(q)
                        attr[a][b]=attr[a][b]+1
print 'success!'
for i in range(len(node)):
    for j in range(len(go)):
        file3.write(str(attr[i][j]))
        file3.write(' ')
    file3.write('\n')
file3.close()

