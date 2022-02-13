import numpy as np
import pandas as pd

def apply_bit_mutation(gene, mutationrate=0.1, randseed=100):
    """
    Apply mutation to a bit gene
    :param gene: 1D numpy array of bit gene
    :param mutationrate: float of mutation rate (0 to 1)
    :param randseed: int number for random state
    :return: mutant gene
    """
    np.random.seed(randseed)
    return gene * (np.random.random(size=len(gene)) <= mutationrate)


def apply_bit_crossover(gene_a, gene_b):
    """
    Cross over bit genes
    :param gene_a: 1d numpy array of bit gene
    :param gene_b: 1d numpy array of bit gene
    :return: 2x 1d numpy array of children genes
    """
    # get genesize
    genesize = len(gene_a)
    # get cut size
    if genesize % 2 == 1:
        cutsize = int((genesize - 1) / 2)
    else:
        cutsize = int(genesize / 2)
    # make copies
    children_a = gene_a.copy()
    children_b = gene_b.copy()
    # cross over copies
    children_a[:cutsize] = gene_b[cutsize:]
    children_b[:cutsize] = gene_a[cutsize:]
    return children_a, children_b


def rws(pscores, samplesize, randseed=100):
    """
    Selection of ids unsing the Roulette Wheel Selection
    :param pscores: 1d numpy array of pscores (must sum up to one)
    :param samplesize: int number of samples
    :param randseed: int number for random state
    :return: 1d numpy array of selected ids
    """
    # generate random arrays
    np.random.seed(randseed)
    _rands = np.random.random(size=samplesize)
    # deploy ids arrays
    _ids = np.zeros(samplesize, dtype='int32')
    # take samples using the wheel
    for i in range(samplesize):
        # reset the accumulated p
        _acc_p = 0
        # reset the id counter
        j = 0
        # run wheel
        while True:
            # accumulate p
            _acc_p = _acc_p + pscores[j]
            # stop condition
            if _acc_p > _rands[i]:
                _ids[i] = j  # set id
                break
            else:
                j = j + 1
                if j >= len(pscores):
                    break
    return _ids


def fuzzy_transition(array, a, b, ascending=True, kind='linear'):
    """

    Fuzzify a numpy array by a transition from a to b values

    :param array: numpy array
    :param a: float initial threshold
    :param b: float terminal threshold
    :param ascending: boolean
    :param kind: string type of fuzzy. options: 'senoid' and 'linear' (trapezoid)
    :return: numpy array
    """
    if a == b:
        transition = (array * 0.0) + 0.5
    else:
        if ascending:
            if kind == 'senoid':
                transition = (array >= b) + (array > a) * (array < b) * (-0.5 * np.cos(np.pi * (array - a)/(b - a)) + 0.5 )
            if kind == 'linear':
                transition = (array >= b) + (array > a) * (array < b) * (( array / (b - a)) - (a / (b - a)))
        else:
            if kind == 'senoid':
                transition = (array <= a) + (array > a) * (array < b) * (0.5 * np.cos(np.pi * (array - a)/(b - a)) + 0.5)
            if kind == 'linear':
                transition = (array <= a) + (array > a) * (array < b) * ((- array / (b - a)) + (b / (b - a)))
    return transition


def generate_population(genesize=10, dnasize=10, popsize=10, randseed=100):
    """
    generate a random population of dnas
    :param genesize: int of gene size, numeber of alleles, or bits
    :param dnasize: int number of genes per DNA (individual or solution)
    :param popsize: int number of DNAs in population
    :param randseed: int number for random state
    :return: 3D numpy array of population
    """
    # deploy random state
    np.random.seed(randseed)
    # deploy blank array
    _pop = np.zeros(shape=(popsize, dnasize, genesize), dtype='int16')
    # insert gene by gene for the sake of memory saving
    for i in range(popsize):
        for j in range(dnasize):
            # deploy random gene
            _lcl_rnd = np.random.random(size=genesize)
            # apply boolean operation
            _lcl_bool = _lcl_rnd > 0.5
            # insert in local row:
            _pop[i][j] = _lcl_bool
    return _pop


def evolve_2d_function(lo, hi,
                       generations=10,
                       popsize=10,
                       genesize=10,
                       mutrate=0.1,
                       elitism=False,
                       trace=False,
                       kind='rastrigin'):
    """

    :param lo:
    :param hi:
    :param generations:
    :param popsize:
    :param genesize:
    :param mutrate:
    :param elitism:
    :param trace:
    :return:
    """
    from datetime import datetime
    from code import decode_binary_float
    if kind == 'rastrigin':
        from benchmark import rastrigin_2d
    elif kind == 'sphere':
        from benchmark import sphere_2d
    elif kind == 'himmelblaus':
        from benchmark import himmelblaus
    else:
        from benchmark import sphere_2d
    # get 2d center
    mid = (hi - lo) / 2

    # set random state
    seed = int(str(datetime.now())[-6:])

    # get seeds
    seeds = np.random.randint(100, 1000, size=generations + 1)

    # declare curve dataframe
    evolution_curve = pd.DataFrame({'Gen': np.arange(0, generations)})
    evolution_curve['Best_S'] = 0.0

    # generate initial population with popsize
    dnasize = 2
    parents = generate_population(genesize=genesize, dnasize=dnasize, popsize=popsize, randseed=seeds[0])
    offspring = parents.copy()  # copy array to use later
    if trace:
        traced = np.zeros(shape=(generations, popsize, dnasize, genesize), dtype='int16')
        traced_scores = np.zeros(shape=(generations, popsize))

    # evaluate initial population fitness scores

    # -- deploy scores array
    scores = np.zeros(popsize, dtype='float16')

    # -- evaluate each solution
    for i in range(popsize):
        #
        # -- decode y and x
        lcl_x = decode_binary_float(gene=parents[i][0], lo=lo, hi=hi)
        lcl_y = decode_binary_float(gene=parents[i][1], lo=lo, hi=hi)
        #
        # -- compute local score
        if kind == 'rastrigin':
            scores[i] = rastrigin_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
        elif kind == 'sphere':
            scores[i] = sphere_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
        elif kind == 'himmelblaus':
            scores[i] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
        else:
            scores[i] = sphere_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)


    # -- evolve generations
    for g in range(generations):
        #print('\n\nGeneration {}'.format(g))
        # store best score
        evolution_curve['Best_S'].values[g] = np.max(scores)
        if trace:
            traced[g] = parents.copy()
            traced_scores[g] = scores.copy()

        # generate offspring using mating pool ids and parents
        for i in range(0, popsize - 1, 2):
            for j in range(dnasize):
                # retrieve genes
                gene_a = parents[i][j]
                gene_b = parents[i + 1][j]
                # insert in mating pool array:
                parents[i][j] = gene_a
                parents[i + 1][j] = gene_b
                # recombine new genes
                new_gene_a, new_gene_b = apply_bit_crossover(gene_a=gene_a, gene_b=gene_b)
                # mutation
                mutant_gene_a = apply_bit_mutation(gene=new_gene_a, mutationrate=mutrate, randseed=seeds[g + 1])
                mutant_gene_b = apply_bit_mutation(gene=new_gene_b, mutationrate=mutrate, randseed=seeds[g + 1])
                # insert in offspring array:
                offspring[i][j] = mutant_gene_a
                offspring[i + 1][j] = mutant_gene_b

        # evaluate offspring
        # -- deploy scores array
        off_scores = np.zeros(popsize, dtype='float16')

        # -- evaluate each solution
        for i in range(popsize):
            #
            # -- decode y and x
            lcl_x = decode_binary_float(gene=offspring[i][0], lo=lo, hi=hi)
            lcl_y = decode_binary_float(gene=offspring[i][1], lo=lo, hi=hi)
            #
            # -- compute local score
            if kind == 'rastrigin':
                off_scores[i] = rastrigin_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
            elif kind == 'sphere':
                off_scores[i] = sphere_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
            elif kind == 'himmelblaus':
                off_scores[i] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
            else:
                off_scores[i] = sphere_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)

        # merge offspring with parents
        pool = np.append(parents, offspring, axis=0)
        pool_scores = np.append(scores, off_scores)

        # select new parents generation from
        # -- selection method
        if elitism:
            # sort pool
            lcl_df = pd.DataFrame({'S': pool_scores})
            lcl_df = lcl_df.sort_values(by='S', ascending=False)
            # slice a popsize sample
            pool_ids = np.array(lcl_df.index)[:popsize]
        else:  # RWS
            # -- compute pscores
            fscores = fuzzy_transition(array=pool_scores,
                                       a=np.min(pool_scores),
                                       b=np.max(pool_scores),
                                       kind='linear',
                                       ascending=True)
            pool_pscores = fscores / np.sum(fscores)
            pool_ids = rws(pscores=pool_pscores, samplesize=popsize, randseed=seeds[g + 1])

        # retrieve parents and scores from pool:
        for i in range(popsize):
            parents[i] = pool[pool_ids[i]]
            scores[i] = pool_scores[pool_ids[i]]
        #print(g)
        #print(scores)
    # return best solution
    last_best_score = np.max(scores)
    for i in range(popsize):
        if scores[i] == last_best_score:
            last_best_score_id = i
    best_x = decode_binary_float(gene=parents[last_best_score_id][0], lo=lo, hi=hi)
    best_y = decode_binary_float(gene=parents[last_best_score_id][1], lo=lo, hi=hi)

    # return data
    out_dct = {'Curve': evolution_curve, 'X': best_x, 'Y': best_y}
    if trace:
        traced_solutions = np.zeros(shape=(generations, dnasize, popsize))
        aux_dct = dict()
        for g in range(generations):
            # express solutions
            for i in range(dnasize):
                for j in range(popsize):
                    traced_solutions[g][i][j] = decode_binary_float(gene=traced[g][j][i], lo=lo, hi=hi)
            aux_dct['G'+ str(g) + '_x'] = traced_solutions[g][0]
            aux_dct['G' + str(g) + '_y'] = traced_solutions[g][1]
            aux_dct['G' + str(g) + '_S'] = traced_scores[g]
        out_dct['Traced'] = pd.DataFrame(aux_dct)

    return out_dct


