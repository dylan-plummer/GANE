file2=open("krogan2006core_node.txt")
file3=open("krogan2006core_go_information.txt","w")


for i in file2:
    i=i.rstrip('\n')
    print i
    file3.write(i)
    file3.write(' ')
    file1 = open("dataset/go_slim_mapping.tab.txt")
    for j in file1:
        node_name = j.split('	')[0]
        go_tag=j.split('	')[3]
        go = j.split('	')[5]
        if node_name==i:
            if go_tag=="P" or go_tag=='F' :
                if go!='' :
                    file3.write(go)
                    file3.write(' ')
    file1.close()
    file3.write('\n')

file2.close()
file3.close()

