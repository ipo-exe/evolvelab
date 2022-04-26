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


def evolve_2d_function(lo_x, hi_x, lo_y, hi_y, mid_x, mid_y,
                       generations=10,
                       popsize=10,
                       genesize=10,
                       mutrate=0.1,
                       elitism=False,
                       trace=False,
                       kind='rastrigin',
                       tui=False):
    """

    2-D evolutionary algorithm.

    Pseudo-code:
    start
        get popsize
        get generations
        get elitism
        generate parents with popsize
        evaluate parents fitness scores
        set g = 0
        repeat until g > generations:
            generate offspring using parents
            evaluate offspring fitness scores
            merge offspring with parents in pool
            select next generation parents from pool
            if elitism == True:
                sort pool by fitness scores
                select the best popsize sample
            else:
                select from pool a popsize sample using RWS method
            set parents = sample
            set g = g + 1
        return parents
    end

    :param lo_x: float of x lowest search space
    :param hi_x: float of x highest search space
    :param lo_y: float of y lowest search space
    :param hi_y: float of y highest search space
    :param mid_x: float center of x
    :param mid_y: float center of y
    :param generations: int number of generations
    :param popsize: int number of population size
    :param genesize: int number of gene bits
    :param mutrate: float 0 to 1 - mutation rate probability
    :param elitism: boolean to use deterministic elitist selection
    :param trace: boolean to trace back the full evolution
    :param kind: string of 2d benchmark function type
    :param tui: boolean to screen display
    :return: output dictionary
    """
    from datetime import datetime
    from backend import status
    from code import decode_binary_float
    if kind == 'rastrigin':
        from benchmark import rastrigin_2d
    elif kind == 'paraboloid':
        from benchmark import paraboloid_2d
    elif kind == 'himmelblaus':
        from benchmark import himmelblaus
    else:
        from benchmark import paraboloid_2d
    # set random state
    seed = int(str(datetime.now())[-6:])

    # get seeds
    seeds = np.random.randint(100, 1000, size=generations + 1)

    # declare curve dataframe
    evolution_curve = pd.DataFrame({'Gen': np.arange(0, generations)})
    evolution_curve['Best_S'] = 0.0
    evolution_curve['p95'] = 0.0
    evolution_curve['p50'] = 0.0
    evolution_curve['p05'] = 0.0

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

        # -- decode y and x
        lcl_x = decode_binary_float(gene=parents[i][0], lo=lo_x, hi=hi_x)
        lcl_y = decode_binary_float(gene=parents[i][1], lo=lo_y, hi=hi_y)

        # -- compute local score
        if kind == 'rastrigin':
            scores[i] = rastrigin_2d(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
        elif kind == 'paraboloid':
            scores[i] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
        elif kind == 'himmelblaus':
            scores[i] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
        else:
            scores[i] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)

    # -- evolve generations
    for g in range(generations):
        if tui:
            status('Generation {}'.format(g))
        # store best score
        evolution_curve['Best_S'].values[g] = np.max(scores)
        evolution_curve['p95'].values[g] = np.percentile(scores, 95)
        evolution_curve['p50'].values[g] = np.percentile(scores, 50)
        evolution_curve['p05'].values[g] = np.percentile(scores, 5)
        if trace:
            traced[g] = parents.copy()
            traced_scores[g] = scores.copy()

        # generate offspring using parents
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
            lcl_x = decode_binary_float(gene=offspring[i][0], lo=lo_x, hi=hi_x)
            lcl_y = decode_binary_float(gene=offspring[i][1], lo=lo_y, hi=hi_y)
            #
            # -- compute local score
            if kind == 'rastrigin':
                off_scores[i] = rastrigin_2d(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
            elif kind == 'paraboloid':
                off_scores[i] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
            elif kind == 'himmelblaus':
                off_scores[i] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid_x, y0=mid_y)
            else:
                off_scores[i] = paraboloid_2d(x=lcl_x, y=lcl_y, xx0=mid_x, y0=mid_y)

        # merge offspring with parents
        pool = np.append(parents, offspring, axis=0)
        pool_scores = np.append(scores, off_scores)

        # select new parents generation from
        # -- selection method
        if elitism:
            # sort pool
            lcl_df = pd.DataFrame({'S': pool_scores})
            lcl_df = lcl_df.sort_values(by='S', ascending=False)
            # slice off a top-popsize sample
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

    # return best solution
    last_best_score = np.max(scores)
    for i in range(popsize):
        if scores[i] == last_best_score:
            last_best_score_id = i
    best_x = decode_binary_float(gene=parents[last_best_score_id][0], lo=lo_x, hi=hi_x)
    best_y = decode_binary_float(gene=parents[last_best_score_id][1], lo=lo_y, hi=hi_y)

    # return data
    out_dct = {'Curve': evolution_curve, 'X': best_x, 'Y': best_y}
    if trace:
        traced_solutions = np.zeros(shape=(generations, dnasize, popsize))
        aux_dct = dict()
        los = [lo_x, lo_y]
        his = [hi_x, hi_y]
        for g in range(generations):
            # express solutions
            for i in range(dnasize):
                for j in range(popsize):
                    traced_solutions[g][i][j] = decode_binary_float(gene=traced[g][j][i], lo=los[i], hi=his[i])
            aux_dct['G'+ str(g) + '_x'] = traced_solutions[g][0]
            aux_dct['G' + str(g) + '_y'] = traced_solutions[g][1]
            aux_dct['G' + str(g) + '_S'] = traced_scores[g]
        out_dct['Traced'] = pd.DataFrame(aux_dct)

    return out_dct


def evolve(genes_df, generations=10, popsize=10, mutt=1, coarse=True, trace=True):

    from benchmark import rastrigin_2d, paraboloid_2d, himmelblaus

    def express_gene(x, y_lo, y_hi, x_hi=255):
        """
        express a gene in float values
        :param x: gene array or int
        :param y_lo: lower bound array or float of y
        :param y_hi: upper bound array or float of y
        :param x_hi: upper bound array or float of x
        :return: float array or float of gene phenotype
        """
        _m = (y_hi - y_lo) / x_hi
        return (_m * x) + y_lo


    def apply_gene_average(gene_a, gene_b, origin_type):
        gene_a = np.array(gene_a, dtype='float32')
        gene_b = np.array(gene_b, dtype='float32')
        output = (gene_a + gene_b) / 2
        output = output.astype(origin_type)
        return output


    from datetime import datetime
    from backend import status
    from sys import getsizeof
    if coarse:
        _dtype = 'uint8'
        _high = 255
    else:
        _dtype = 'uint16'
        _high = 65535

    # get gene size:
    genesize = len(genes_df['Hi'])

    # set overall random state
    seed = int(str(datetime.now())[-6:])
    # get random seeds for each generation
    seeds = np.random.randint(100, 1000, size=generations)

    # declare curve dataframe
    evolution_curve = pd.DataFrame({'Gen': np.arange(0, generations)})
    evolution_curve['Best_S'] = 0.0
    evolution_curve['p95'] = 0.0
    evolution_curve['p50'] = 0.0
    evolution_curve['p05'] = 0.0

    # set tracing variables
    if trace:
        traced_parents = np.zeros(shape=(generations, popsize, genesize), dtype=_dtype)
        traced_parents_scores = np.zeros(shape=(generations, popsize), dtype='float16')
        traced_offpring = np.zeros(shape=(generations, popsize, genesize), dtype=_dtype)
        traced_offpring_scores = np.zeros(shape=(generations, popsize), dtype='float16')
        status(msg='traced parents array size: {} KB'.format(getsizeof(traced_parents) / 1000), process=False)
        status(msg='traced scores array size: {} KB'.format(getsizeof(traced_parents_scores) / 1000), process=False)
        status(msg='traced parents array size: {} KB'.format(getsizeof(traced_offpring) / 1000), process=False)
        status(msg='traced scores array size: {} KB'.format(getsizeof(traced_offpring_scores) / 1000), process=False)

    # get initial population:
    parents = np.random.randint(low=0, high=_high, size=(popsize, genesize), dtype=_dtype)
    parents = int(_high / 2) + np.zeros(shape=(popsize, genesize), dtype=_dtype)

    # generations loop:
    for g in range(0, generations):
        status('generation {}'.format(g))

        # reset random state
        np.random.seed(seeds[g])

        # shuffle parents
        np.random.shuffle(parents)

        # get offspring from parents:
        offspring = 0 * parents.copy()
        for i in range(len(parents)):
            if i < len(parents) - 1:
                _first_id = i
                _secnd_id = i + 1
            else:
                _first_id = i
                _secnd_id = 0
            # reproduce genes by gene average
            lcl_gene_offs = apply_gene_average(gene_a=parents[_first_id], gene_b=parents[_secnd_id], origin_type=_dtype)
            offspring[i] = lcl_gene_offs.astype(_dtype)

        # apply variation operator:
        mutation_mask = np.random.normal(loc=0, scale=mutt, size=np.shape(offspring)).astype(dtype=_dtype)
        offspring = offspring + mutation_mask

        # evaluate full population
        parents_scores = np.zeros(shape=popsize, dtype='float')
        offspring_scores = np.zeros(shape=popsize, dtype='float')
        # declare a map (dictionary) to assess the full population
        population = {'Parents': {'Genes': parents, 'Scores': parents_scores},
                      'Offspring': {'Genes': offspring, 'Scores': offspring_scores}}
        for i in range(2 * popsize):
            # decide key
            lcl_key = 'Parents'
            lcl_id = i
            if i >= popsize:
                lcl_key = 'Offspring'
                lcl_id = i - popsize
            # express parent gene
            lcl_gene = express_gene(x=population[lcl_key]['Genes'][lcl_id],
                                    x_hi=_high,
                                    y_lo=genes_df['Lo'].values,
                                    y_hi=genes_df['Hi'].values)
            #
            #
            # compute objective function
            ## population[lcl_key]['Scores'][lcl_id] = rastrigin_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            population[lcl_key]['Scores'][lcl_id] = himmelblaus(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            #
            #
            #
        # retrieve best next parents
        scores = np.concatenate((parents_scores, offspring_scores))  # merge scores
        scores_ids = np.arange(len(scores)) # create ids
        scores_df = pd.DataFrame({'Id': scores_ids, 'Score': scores})
        scores_df['New'] = 0.0
        scores_df['New'] = (scores_df['Score'].values * (scores_df['Score'].values < 90)) + (100 * (scores_df['Score'].values >= 90))
        # todo here
        scores_df.sort_values(by='New', ascending=False, inplace=True)  # sort scores and ids

        for i in range(popsize):
            # decide key
            lcl_score_id = scores_df['Id'].values[i]
            if lcl_score_id >= popsize:
                lcl_key = 'Offspring'
                lcl_id = lcl_score_id - popsize
            else:
                lcl_key = 'Parents'
                lcl_id = lcl_score_id
            #
            # replace parent
            parents[i] = population[lcl_key]['Genes'][lcl_id]

        # store best score
        evolution_curve['Best_S'].values[g] = np.max(parents_scores)
        evolution_curve['p95'].values[g] = np.percentile(parents_scores, 95)
        evolution_curve['p50'].values[g] = np.percentile(parents_scores, 50)
        evolution_curve['p05'].values[g] = np.percentile(parents_scores, 5)

        # trace
        if trace:
            traced_parents[g] = parents.copy()
            traced_parents_scores[g] = parents_scores.copy()
            traced_offpring[g] = offspring.copy()
            traced_offpring_scores[g] = offspring_scores.copy()

    last_best_solution = express_gene(x=parents[0],
                                      x_hi=_high,
                                      y_lo=genes_df['Lo'].values,
                                      y_hi=genes_df['Hi'].values)
    out_dct = {'Curve': evolution_curve}
    for i in range(len(genes_df['Labels'].values)):
        out_dct['Best_{}'.format(genes_df['Labels'].values[i])] = last_best_solution[i]

    if trace:
        aux_generations = np.zeros(generations * popsize, dtype='uint16')
        trace_df = pd.DataFrame({'Gen': aux_generations})
        trace_df['Score'] = 0.0
        for k in genes_df['Labels'].values:
            trace_df[k] = 0.0
        counter = 0
        for g in range(generations):
            for i in range(popsize):
                lcl_solution = express_gene(x=traced_parents[g][i],
                                            x_hi=_high,
                                            y_lo=genes_df['Lo'].values,
                                            y_hi=genes_df['Hi'].values)
                trace_df['Gen'].values[counter] = g
                trace_df['Score'].values[counter] = traced_parents_scores[g][i]
                for j in range(len(genes_df['Labels'].values)):
                    trace_df[genes_df['Labels'].values[j]].values[counter] = lcl_solution[j]
                counter = counter + 1
        out_dct['Traced'] = trace_df
    return out_dct

'''
ranges_df = pd.DataFrame({'Lo': [-20, -20],
                          'Hi': [2, 2],
                          'Labels': ['X', 'Y']})

e = evolve(genes_df=ranges_df,
           popsize=5,
           generations=3,
           coarse=False,
           trace=True,
           mutt=1)
print(e)

import matplotlib.pyplot as plt

curve = e['Curve']

plt.plot(curve['Gen'], curve['Best_S'])
plt.plot(curve['Gen'], curve['p95'])
plt.plot(curve['Gen'], curve['p50'])
plt.plot(curve['Gen'], curve['p05'])
plt.ylim(0, 101)
plt.show()


'''


