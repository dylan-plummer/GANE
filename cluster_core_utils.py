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