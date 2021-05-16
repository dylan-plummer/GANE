import os
from numpy import *

from numpy.linalg import *

def f_key(a):
    return(a[-1])


def density_score(temp_set, matrix):
    temp_density_score = 0.
    for m in temp_set:
        for n in temp_set:
            if n != m and matrix[m, n] != 0:
                temp_density_score += matrix[m, n]

    temp_density_score = temp_density_score / (len(temp_set) * (len(temp_set) - 1))
    return temp_density_score


def merge_cliques(new_cliques_set, matrix):
    seed_clique = []

    while (True):
        temp_cliques_set = []
        if len(new_cliques_set) >= 2:
            seed_clique.append(new_cliques_set[0])

            for i in range(1, len(new_cliques_set)):
                if len(new_cliques_set[i].intersection(new_cliques_set[0])) == 0:
                    temp_cliques_set.append(new_cliques_set[i])
                elif len(new_cliques_set[i].difference(new_cliques_set[0])) >= 3:
                    temp_cliques_set.append(new_cliques_set[i].difference(new_cliques_set[0]))

            cliques_set = []

            for i in temp_cliques_set:

                clique_score = density_score(i, matrix)
                temp_list = []
                for j in i:
                    temp_list.append(j)

                temp_list.append(clique_score)
                cliques_set.append(temp_list)

            cliques_set.sort(key=f_key, reverse=True)

            new_cliques_set = []
            for i in range(len(cliques_set)):
                temp_set = set([])
                for j in range(len(cliques_set[i]) - 1):
                    temp_set.add(cliques_set[i][j])
                new_cliques_set.append(temp_set)

        elif len(new_cliques_set) == 1:
            seed_clique.append(new_cliques_set[0])
            break
        else:
            break

    return seed_clique


def expand_cluster(seed_clique, all_protein_set, matrix, expand_thres):
    expand_set = []
    complex_set = []

    for instance in seed_clique:
        avg_node_score = density_score(instance, matrix)

        temp_set = set([])
        for j in all_protein_set.difference(instance):
            temp_score = 0.
            for n in instance:
                temp_score += matrix[n, j]
            temp_score /= len(instance)

            if (temp_score ) >= expand_thres:
                temp_set.add(j)
        expand_set.append(temp_set)
    for i in range(len(seed_clique)):
        complex_set.append(seed_clique[i].union(expand_set[i]))

    return (complex_set)


def cluster(dataset, dset_file, method, protein_out_file='protein.temp', cliques_file="cliques", ppi_pair_file="ppi.pair", ppi_matrix_file="ppi.matrix"):
    f = open("plots/%s/%s/%s_sim.txt" % (dataset, method, dataset), "r")
    f_protein_out = open(protein_out_file, "w")
    Dic_map = {}
    index = 0
    Node1 = []
    Node2 = []
    Weight = []
    All_node = set([])
    All_node_index = set([])

    for line in f:
        line = line.strip().split()
        if len(line) == 3:
            Node1.append(line[0])
            All_node.add(line[0])
            Node2.append(line[1])
            All_node.add(line[1])
            Weight.append(float(line[2]))
            if line[0] not in Dic_map:
                Dic_map[line[0]] = index
                f_protein_out.write(line[0] + "\n")
                All_node_index.add(index)
                index += 1
            if line[1] not in Dic_map:
                Dic_map[line[1]] = index
                f_protein_out.write(line[1] + "\n")
                All_node_index.add(index)
                index += 1
    Node_count = index
    f.close()
    f_protein_out.close()
    # print Dic_map

    f.close()

    ######dic_map to map_dic###########
    Map_dic = {}
    for key in list(Dic_map.keys()):
        Map_dic[Dic_map[key]] = key

    # print Map_dic

    ######bulid Adj_matrix###########

    Adj_Matrix = mat(zeros((Node_count, Node_count), dtype=float))

    if len(Node1) == len(Node2):

        for i in range(len(Node1)):
            if Node1[i] in Dic_map and Node2[i] in Dic_map:
                Adj_Matrix[Dic_map[Node1[i]], Dic_map[Node2[i]]] = Weight[i]
                Adj_Matrix[Dic_map[Node2[i]], Dic_map[Node1[i]]] = Weight[i]
    # print Adj_Matrix.shape[0]

    os.system(
        "ConvertPPI.exe " + "dataset/%s.txt" % dset_file + " " + protein_out_file + " " + ppi_pair_file + " " + ppi_matrix_file)
    os.system(
        "Mining_Cliques.exe " + ppi_matrix_file + " " + "1" + " " + "3" + " " + str(Node_count) + " " + cliques_file)
    cliques_set = []
    f = open(cliques_file, "r")
    for line in f:
        temp_set = []
        line = line.strip().split()
        for i in range(1, len(line)):
            temp_set.append(int(line[i]))
        cliques_set.append(temp_set)

    f.close()
    avg_clique_score = 0.

    for instance in cliques_set:
        clique_score = density_score(instance, Adj_Matrix)
        avg_clique_score += clique_score
        instance.append(clique_score)
    avg_clique_score /= len(cliques_set)
    cliques_set.sort(key=f_key, reverse=True)

    new_cliques_set = []
    for i in range(len(cliques_set)):
        temp_set = set([])
        for j in range(len(cliques_set[i]) - 1):
            temp_set.add(cliques_set[i][j])
        new_cliques_set.append(temp_set)
    #print('merging cliques...')
    seed_clique = merge_cliques(new_cliques_set, Adj_Matrix)
    #print('expanding clusters...')
    expand_thres = 0.3
    complex_set = expand_cluster(seed_clique, All_node_index, Adj_Matrix, expand_thres)
    #print("##########output predicted complexes##########\n")
    final_file = open("result/final_%s_attr_output" % dataset, "w")

    for i in range(len(complex_set)):

        line = ""
        for m in complex_set[i]:
            line += Map_dic[m] + " "
        line += "\n"

        final_file.write(line)
    final_file.close()


def eval(dataset):
    str = "result/final_%s_attr_output" % dataset
    file = open(str)
    file1 = open("dataset/Form_CYC20083.txt")
    #file1=open("dataset/golden_standard.txt")

    # g=nx.Graph()
    #predicted_num = len(file.readlines())
    #reference_num = len(file1.readlines())
    reference_complex = []
    for j in file1:
        j = j.rstrip()
        j = j.rstrip('\n')
        complex_list = j.split('	')
        reference_complex.append(complex_list)

    predicted_complex = []
    for i in file:
        i = i.rstrip()
        i = i.rstrip('\n')
        node_list = i.split(' ')
        predicted_complex.append(node_list)

    predicted_num = len(predicted_complex)
    reference_num = len(reference_complex)
    # precision
    number = 0
    c_number = 0
    row = 1
    for i in predicted_complex:
        overlapscore = 0.0
        for j in reference_complex:
            set1 = set(i)
            set2 = set(j)
            overlap = set1 & set2
            score = float((pow(len(overlap), 2))) / float((len(set1) * len(set2)))
            if (score > overlapscore):
                overlapscore = score
        if (overlapscore > 0.2):
            number = number + 1
            #print(row, end=' ')
            #print(" ", end=' ')
        row = row + 1
    # recall
    for i in reference_complex:
        overlapscore = 0.0
        for j in predicted_complex:
            set1 = set(i)
            set2 = set(j)
            overlap = set1 & set2
            score = float((pow(len(overlap), 2))) / float((len(set1) * len(set2)))
            if (score > overlapscore):
                overlapscore = score
        if (overlapscore > 0.25):
            c_number = c_number + 1
    # sn
    T_sum1 = 0.0
    N_sum = 0.0
    for i in reference_complex:
        max = 0.0
        for j in predicted_complex:
            set1 = set(i)
            set2 = set(j)
            overlap = set1 & set2
            if len(overlap) > max:
                max = len(overlap)
        T_sum1 = T_sum1 + max
        N_sum = N_sum + len(set1)
    # ppv
    T_sum2 = 0.0
    T_sum = 0.0
    for i in predicted_complex:
        max = 0.0
        for j in reference_complex:
            set1 = set(i)
            set2 = set(j)
            overlap = set1 & set2
            T_sum = T_sum + len(overlap)
            if len(overlap) > max:
                max = len(overlap)
        T_sum2 = T_sum2 + max

    #print("\n")
    print(number, predicted_num)  # matched predicted complex number
    # print c_number,reference_num# matched reference complex number
    precision = float(number / float(predicted_num))
    recall = float(c_number / float(reference_num))
    F1 = float((2 * precision * recall) / (precision + recall))
    Sn = float(T_sum1) / float(N_sum)
    PPV = float(T_sum2) / float(T_sum)
    Acc = pow(float(Sn * PPV), 0.5)
    print("precision", end=' ')
    print(precision)
    print("recall", end=' ')
    print(recall)
    print("F1", end=' ')
    print(F1)
    print("Sn", end=' ')
    print(Sn)
    print("PPV", end=' ')
    print(PPV)
    print("Acc", end=' ')
    print(Acc)
    file.close()
    file1.close()
    return precision, recall, F1, Sn, PPV, Acc