import numpy as np

str1="dataset/biogrid.txt"
str2="dataset/biogrid_attr_vector224.txt"
str3="dataset/biogrid_attr_sim224.txt"

file=open("H32")
file1=open(str1)
file2=open(str2,'w')
print "get the vector representation: "
node=[]
for j in file1:
    temp1=j.split('	')[0]
    temp2=j.split('	')[1].rstrip('\n')
    if temp1 not in node:
        node.append(temp1)
    if temp2 not in node:
        node.append(temp2)
file1.close()

d=[]
for i in file:
    d.append(i)

for i in range(len(node)):
    file2.write(node[i])
    file2.write(' ')
    file2.write(d[i])
    file2.write('\n')
file2.close()
print "calculate the similarity between two nodes:"
file=open(str2)
file1=open(str1)
file2=open(str3,'w')

def cos_sim(vector1,vector2):
    dot_product = 0.0
    normA = 0.0
    normB = 0.0
    for a, b in zip(vector1, vector2):
        dot_product += a * b
        normA += a ** 2
        normB += b ** 2
    result=dot_product / ((normA * normB) ** 0.5)
    return result

edge_name_name=[]
for i in file1:
    node_name1=i.split('	')[0]
    node_name2=i.split('	')[1]
    node_name2 = node_name2.split('\n')[0]
    d={}
    d['node_name1']=node_name1
    d['node_name2']=node_name2
    edge_name_name.append(d)

vector=[]
for i in file:
    if not i.strip(): continue
    node_name=i.split(' ',1)[0]
    node_vector=i.split(' ',1)[1].rstrip('\n')
    node_vector = node_vector.split('	')
    node_vector = map(float, node_vector)
    d = {}
    d['node_name'] = node_name
    d['node_vector'] = node_vector
    vector.append(d)

v1=[]
v2=[]
for i in edge_name_name:
    temp1=0
    temp2=0
    for j in vector:
        if(i['node_name1']==j['node_name']):
            v1=np.array(j['node_vector'])
            temp1=1
    for z in vector:
        if(i['node_name2']==z['node_name']):
            v2=np.array(z['node_vector'])
            temp2=1
    if(temp1==1)and(temp2==1):
        result=cos_sim(v1,v2)
        file2.write(i['node_name1'])
        file2.write(' ')
        file2.write(i['node_name2'])
        file2.write(' ')
        file2.write(str(result))
        file2.write('\n')

file.close()
file1.close()
file2.close()