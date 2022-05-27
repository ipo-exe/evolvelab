import numpy as np
import pandas as pd


def evolve(df_dvars,
           n_generations=10,
           n_popsize=10,
           r_std=0.01,
           kind = 'paraboloid',
           b_coarse=True,
           b_trace=True,
           b_explore=False,
           b_minimize=False,
           upper=90,
           lower=80):
    """

    evolutionary algorithm

    :param df_dvars: pandas dataframe of decision variables, with fields: 'Lo', 'Hi' and 'Labels'
    :param n_generations: integer number of generations
    :param n_popsize: integer number of population
    :param r_std: float ratio of standard deviation
    :param kind: string code for objetive function type.
    Options: 'paraboloid', 'rastrigin', 'himmelblaus', 'mixed', 'custom'
    :param b_coarse: boolean for grid resolution definition. True = 0 to 255. False = 0 to 65535
    :param b_trace:  boolean for tracing back all evolution
    :param b_explore: boolean to explore instead of fitting
    :param b_minimize: boolean to minimize instead of maximize
    :param upper: float - upper exploration value of objective function
    :param lower: float - lower exploration value of objective function
    :return: evolution dictionary
    """

    from benchmark import rastrigin_2d, paraboloid_2d, himmelblaus, griewank_2d
    from datetime import datetime
    from backend import status
    from sys import getsizeof


    def express_gene(gene, feno_lo, feno_hi, gene_hi=255):
        """
        express a gene in float values
        :param gene: gene array or int
        :param feno_lo: lower bound array or float of y
        :param feno_hi: upper bound array or float of y
        :param gene_hi: upper bound array or float of x
        :return: float array or float of gene phenotype
        """
        _m = (feno_hi - feno_lo) / gene_hi
        return (_m * gene) + feno_lo


    def compute_fitness(kind, lcl_gene):
        if kind == 'rastrigin':
            lcl_score_value = rastrigin_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
        elif kind == 'himmelblaus':
            lcl_score_value = himmelblaus(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
        elif kind == 'griewank':
            lcl_score_value = griewank_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
        elif kind == 'paraboloid':
            lcl_score_value = paraboloid_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
        elif kind == 'mixed':
            # rast
            lcl_score_value1 = rastrigin_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
            # himm
            lcl_score_value2 = himmelblaus(x=lcl_gene[2], y=lcl_gene[3], x0=0, y0=0, level=100)
            # grie
            lcl_score_value3 = griewank_2d(x=lcl_gene[4], y=lcl_gene[5], x0=0, y0=0, level=100)
            # parab
            lcl_score_value4 = paraboloid_2d(x=lcl_gene[6], y=lcl_gene[7], x0=0, y0=0, level=100)
            lcl_score_value = np.average([lcl_score_value1, lcl_score_value2, lcl_score_value3, lcl_score_value4])
        elif kind == 'custom':
            lcl_score_value = 100 # custom function
        else:
            lcl_score_value = paraboloid_2d(x=lcl_gene[0], y=lcl_gene[1], x0=0, y0=0, level=100)
        return lcl_score_value


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
    n_genesize = len(df_dvars['Hi'])

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
    ##grd_parents = int(n_high / 10) + np.zeros(shape=(n_popsize, n_genesize), dtype=s_dtype)

    # generations loop:
    for g in range(0, n_generations):
        status('generation {}'.format(g))

        # reset random state
        np.random.seed(v_seeds[g])

        # get offspring from parents:
        grd_offspring = grd_parents.copy()

        # recombination idea:
        """
        # shuffle parents
        np.random.shuffle(grd_parents)
        
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
        """

        # apply variation operator (mutation -- rizome propagation):
        grd_muttation_delta = np.random.normal(loc=0, scale=std, size=np.shape(grd_offspring)).astype(dtype=s_dtype)
        grd_offspring = grd_offspring + grd_muttation_delta # offsprint +- a normal change

        # evaluate full population at once

        # declare scores arrays
        if g == 0:
            v_parents_scores = np.zeros(shape=n_popsize, dtype='float') # initialize parents scores array
        v_offspring_scores = np.zeros(shape=n_popsize, dtype='float')
        # declare a map (dictionary) to assess the full population without duplicate variables in memory
        dct_population = {'Parents': {'Genes': grd_parents, 'Scores': v_parents_scores},
                          'Offspring': {'Genes': grd_offspring, 'Scores': v_offspring_scores}}

        # evaluation operator loop
        if g == 0: # evaluate parents and offspring
            for i in range(2 * n_popsize):
                # decide key
                lcl_key = 'Parents'
                lcl_id = i
                if i >= n_popsize:
                    lcl_key = 'Offspring'
                    lcl_id = i - n_popsize
                # express parent gene
                lcl_gene = express_gene(gene=dct_population[lcl_key]['Genes'][lcl_id],
                                        gene_hi=n_high,
                                        feno_lo=df_dvars['Lo'].values,
                                        feno_hi=df_dvars['Hi'].values)
                # compute objective function
                dct_population[lcl_key]['Scores'][lcl_id] = compute_fitness(kind=kind, lcl_gene=lcl_gene)
        else:  # evaluate only offspring
            for i in range(n_popsize):
                # decide key
                lcl_key = 'Offspring'
                lcl_id = i
                # express parent gene
                lcl_gene = express_gene(gene=dct_population[lcl_key]['Genes'][lcl_id],
                                        gene_hi=n_high,
                                        feno_lo=df_dvars['Lo'].values,
                                        feno_hi=df_dvars['Hi'].values)
                # compute objective function
                dct_population[lcl_key]['Scores'][lcl_id] = compute_fitness(kind=kind, lcl_gene=lcl_gene)

        # trace parents and offspring at this point
        if b_trace:
            grd3_tcd_parents[g] = grd_parents.copy()
            grd_tcd_parents_scores[g] = v_parents_scores.copy()
            grd3_tcd_offpring[g] = grd_offspring.copy()
            grd_tcd_offpring_scores[g] = v_offspring_scores.copy()

        # retrieve best next parents
        v_scores = np.concatenate((v_parents_scores, v_offspring_scores))  # merge scores parents FIRST
        df_scores = pd.DataFrame({'Id': np.arange(len(v_scores)), 'Score': v_scores})  # create indexed dataframe

        # exploration or fitness
        if b_explore:
            # Exploration field defined as random integer:
            df_scores['Exploration'] = np.random.randint(0, 100, size=len(df_scores))
            # set Exploration = 0 when not meet the exploration criteria:
            df_scores['Exploration'] = df_scores['Exploration'].values * (df_scores['Score'].values >= lower) * \
                                       ((df_scores['Score'].values <= upper))
            # sort Exploration field:
            df_scores.sort_values(by='Exploration', ascending=b_minimize, inplace=True)  # sort scores and ids
        else: # fitness MAXIMIZING
            # sort Score field:
            df_scores.sort_values(by='Score', ascending=b_minimize, inplace=True)  # sort scores and ids

        # smart collector of ids
        for i in range(n_popsize):

            # use id to retrieve from population
            lcl_score_id = df_scores['Id'].values[i]

            # smart trick to access population
            if lcl_score_id >= n_popsize:  # is an offspring solution
                lcl_key = 'Offspring'
                lcl_id = lcl_score_id - n_popsize
            else:  # is a parent solution
                lcl_key = 'Parents'
                lcl_id = lcl_score_id

            # replace parent:
            grd_parents[i] = dct_population[lcl_key]['Genes'][lcl_id]
            # replace parent score:
            v_parents_scores[i] = df_scores['Score'].values[i]

        # store best score
        df_evolution_curve['Best_S'].values[g] = np.max(v_parents_scores)
        df_evolution_curve['p95'].values[g] = np.percentile(v_parents_scores, 95)
        df_evolution_curve['p50'].values[g] = np.percentile(v_parents_scores, 50)
        df_evolution_curve['p05'].values[g] = np.percentile(v_parents_scores, 5)

    last_best_solution = express_gene(gene=grd_parents[0],
                                      gene_hi=n_high,
                                      feno_lo=df_dvars['Lo'].values,
                                      feno_hi=df_dvars['Hi'].values)
    # define output dict
    dct_output = {'Curve': df_evolution_curve}
    # append keys
    for i in range(len(df_dvars['Labels'].values)):
        dct_output['Best_{}'.format(df_dvars['Labels'].values[i])] = last_best_solution[i]
    # append the evolution tracind
    if b_trace:
        # define dataframe. set Gen field
        df_trace = pd.DataFrame({'Gen': np.zeros(n_generations * n_popsize, dtype='uint16')})
        # set scores field
        df_trace['Score'] = 0.0
        # append more fields
        for k in df_dvars['Labels'].values:
            df_trace[k] = 0.0
        counter = 0
        for g in range(n_generations):
            for i in range(n_popsize):
                # append generation
                df_trace['Gen'].values[counter] = g
                # append score
                df_trace['Score'].values[counter] = grd_tcd_parents_scores[g][i]
                # express local solution
                lcl_solution = express_gene(gene=grd3_tcd_parents[g][i],
                                            gene_hi=n_high,
                                            feno_lo=df_dvars['Lo'].values,
                                            feno_hi=df_dvars['Hi'].values)
                # append solution in fields
                for j in range(len(df_dvars['Labels'].values)):
                    df_trace[df_dvars['Labels'].values[j]].values[counter] = lcl_solution[j]
                counter = counter + 1
        # append to output dictionary
        dct_output['Traced'] = df_trace
    return dct_output

