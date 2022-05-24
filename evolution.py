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


def evolve(df_genes,
           n_generations=10,
           n_popsize=10,
           r_std=0.01,
           r_mutt=0.05,
           kind = 'paraboloid',
           b_coarse=True,
           b_trace=True,
           b_recomb=True,
           b_explore=False,
           upper=90,
           lower=80):

    from benchmark import rastrigin_2d, paraboloid_2d, himmelblaus, griewank_2d
    from datetime import datetime
    from backend import status
    from sys import getsizeof


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


    # resolution setup
    if b_coarse:
        s_dtype = 'uint8'
        n_high = 255
    else:
        s_dtype = 'uint16'
        n_high = 65535
    # standard deviation setup
    std = int(r_std * n_high)

    # get gene size:
    n_genesize = len(df_genes['Hi'])

    # set overall random state
    seed = int(str(datetime.now())[-6:])
    # get random seeds for each generation
    v_seeds = np.random.randint(100, 1000, size=n_generations)

    # declare curve dataframe
    df_evolution_curve = pd.DataFrame({'Gen': np.arange(0, n_generations)})
    df_evolution_curve['Best_S'] = 0.0
    df_evolution_curve['p95'] = 0.0
    df_evolution_curve['p50'] = 0.0
    df_evolution_curve['p05'] = 0.0

    # set tracing variables
    if b_trace:
        grd3_tcd_parents = np.zeros(shape=(n_generations, n_popsize, n_genesize), dtype=s_dtype)
        grd_tcd_parents_scores = np.zeros(shape=(n_generations, n_popsize), dtype='float16')
        grd3_tcd_offpring = np.zeros(shape=(n_generations, n_popsize, n_genesize), dtype=s_dtype)
        grd_tcd_offpring_scores = np.zeros(shape=(n_generations, n_popsize), dtype='float16')
        status(msg='traced parents array size: {} KB'.format(getsizeof(grd3_tcd_parents) / 1000), process=False)
        status(msg='traced scores array size: {} KB'.format(getsizeof(grd_tcd_parents_scores) / 1000), process=False)
        status(msg='traced parents array size: {} KB'.format(getsizeof(grd3_tcd_offpring) / 1000), process=False)
        status(msg='traced scores array size: {} KB'.format(getsizeof(grd_tcd_offpring_scores) / 1000), process=False)

    # get initial population:
    grd_parents = np.random.randint(low=0, high=n_high, size=(n_popsize, n_genesize), dtype=s_dtype)
    ##grd_parents = int(n_high / 2) + np.zeros(shape=(n_popsize, n_genesize), dtype=s_dtype)

    # generations loop:
    for g in range(0, n_generations):
        status('generation {}'.format(g))

        # reset random state
        np.random.seed(v_seeds[g])

        # shuffle parents
        np.random.shuffle(grd_parents)

        # get offspring from parents:
        grd_offspring = grd_parents.copy()

        # apply variation operator (recombination)
        if b_recomb:
            for i in range(len(grd_parents)):
                # smart index selector
                if i < len(grd_parents) - 1:
                    n_1st_id = i
                    n_2nd_id = i + 1
                else:
                    n_1st_id = i
                    n_2nd_id = 0
                # reproduce genes by gene averaging
                lcl_gene_offs = apply_gene_average(gene_a=grd_parents[n_1st_id],
                                                   gene_b=grd_parents[n_2nd_id],
                                                   origin_type=s_dtype)
                # replace in grid
                grd_offspring[i] = lcl_gene_offs.astype(s_dtype)

        # apply variation operator (mutation -- rizome propagation):
        grd_muttation_mask = np.random.normal(loc=0, scale=std, size=np.shape(grd_offspring)).astype(dtype=s_dtype)
        grd_offspring = grd_offspring + grd_muttation_mask # offsprint +- a normal change

        # apply variation operation (mutation -- spore propagation)
        n_spores_determ = int(n_popsize * r_mutt)  # deterministic number of mutations
        n_spores_stocas = int(np.random.normal(loc=n_spores_determ, scale=n_spores_determ / 20)) # normal variation
        grd_new_random_genes = np.random.randint(0, high=n_high, size=n_genesize) # new random genes
        v_muttation_ids = np.random.randint(0, high=n_popsize - 1, size=n_spores_stocas)  # which genes will mutate
        for i in range(n_spores_stocas):
            lcl_id = v_muttation_ids[i]
            grd_offspring[lcl_id] = grd_new_random_genes

        # evaluate full population at once
        v_parents_scores = np.zeros(shape=n_popsize, dtype='float')
        v_offspring_scores = np.zeros(shape=n_popsize, dtype='float')
        # declare a map (dictionary) to assess the full population without duplicate variables in memory
        dct_population = {'Parents': {'Genes': grd_parents, 'Scores': v_parents_scores},
                          'Offspring': {'Genes': grd_offspring, 'Scores': v_offspring_scores}}

        # evaluation operator loop
        for i in range(2 * n_popsize):
            # decide key
            lcl_key = 'Parents'
            lcl_id = i
            if i >= n_popsize:
                lcl_key = 'Offspring'
                lcl_id = i - n_popsize
            # express parent gene
            lcl_gene = express_gene(x=dct_population[lcl_key]['Genes'][lcl_id],
                                    x_hi=n_high,
                                    y_lo=df_genes['Lo'].values,
                                    y_hi=df_genes['Hi'].values)


            # compute objective function
            if kind == 'rastrigin':
                lcl_score_value = rastrigin_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            elif kind == 'himmelblaus':
                lcl_score_value = himmelblaus(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            elif kind == 'griewank':
                lcl_score_value = griewank_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            elif kind == 'paraboloid':
                lcl_score_value = paraboloid_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            else:
                lcl_score_value = paraboloid_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            dct_population[lcl_key]['Scores'][lcl_id] = lcl_score_value
            ##dct_population[lcl_key]['Scores'][lcl_id] =

        # trace parents and offspring at this point
        if b_trace:
            grd3_tcd_parents[g] = grd_parents.copy()
            grd_tcd_parents_scores[g] = v_parents_scores.copy()
            grd3_tcd_offpring[g] = grd_offspring.copy()
            grd_tcd_offpring_scores[g] = v_offspring_scores.copy()

        # retrieve best next parents
        v_scores = np.concatenate((v_parents_scores, v_offspring_scores))  # merge scores parents FIRST
        first_key = 'Parents'
        second_key = 'Offspring'
        df_scores = pd.DataFrame({'Id': np.arange(len(v_scores)), 'Score': v_scores})

        # exploration or fitness
        if b_explore:
            df_scores = pd.DataFrame({'Id': np.arange(len(v_scores)), 'Score': v_scores})
            df_scores['Exploration'] = np.random.randint(0, 100, size=len(df_scores))
            df_scores['Exploration'] = df_scores['Exploration'].values * (df_scores['Score'].values >= lower) * \
                                       ((df_scores['Score'].values <= upper))
            df_scores.sort_values(by='Exploration', ascending=False, inplace=True)  # sort scores and ids
        else: # fitness MAXIMIZING
            # sort scores --
            df_scores.sort_values(by='Score', ascending=False, inplace=True)  # sort scores and ids

        # smart collector of ids
        for i in range(n_popsize):
            # use id to retrieve from population
            lcl_score_id = df_scores['Id'].values[i]
            # smart trick to access population
            if lcl_score_id >= n_popsize:
                lcl_key = second_key
                lcl_id = lcl_score_id - n_popsize
            else:
                lcl_key = first_key
                lcl_id = lcl_score_id
            # replace parent
            grd_parents[i] = dct_population[lcl_key]['Genes'][lcl_id]

        # store best score
        df_evolution_curve['Best_S'].values[g] = np.max(v_parents_scores)
        df_evolution_curve['p95'].values[g] = np.percentile(v_parents_scores, 95)
        df_evolution_curve['p50'].values[g] = np.percentile(v_parents_scores, 50)
        df_evolution_curve['p05'].values[g] = np.percentile(v_parents_scores, 5)



    last_best_solution = express_gene(x=grd_parents[0],
                                      x_hi=n_high,
                                      y_lo=df_genes['Lo'].values,
                                      y_hi=df_genes['Hi'].values)
    # define output dict
    dct_output = {'Curve': df_evolution_curve}
    # append keys
    for i in range(len(df_genes['Labels'].values)):
        dct_output['Best_{}'.format(df_genes['Labels'].values[i])] = last_best_solution[i]
    # append more
    if b_trace:
        # array
        aux_generations = np.zeros(n_generations * n_popsize, dtype='uint16')
        # define dataframe. set Gen field
        df_trace = pd.DataFrame({'Gen': aux_generations})
        # set scores field
        df_trace['Score'] = 0.0
        # append more fields
        for k in df_genes['Labels'].values:
            df_trace[k] = 0.0
        counter = 0
        for g in range(n_generations):
            for i in range(n_popsize):
                # append generation
                df_trace['Gen'].values[counter] = g
                # append score
                df_trace['Score'].values[counter] = grd_tcd_parents_scores[g][i]
                # express local solution
                lcl_solution = express_gene(x=grd3_tcd_parents[g][i],
                                            x_hi=n_high,
                                            y_lo=df_genes['Lo'].values,
                                            y_hi=df_genes['Hi'].values)
                # append solution in fields
                for j in range(len(df_genes['Labels'].values)):
                    df_trace[df_genes['Labels'].values[j]].values[counter] = lcl_solution[j]
                counter = counter + 1
        dct_output['Traced'] = df_trace
    return dct_output

