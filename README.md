# evolvelab
Evolutionary computing lab. For evolutionary computing stuff and experiments.

## The evolutionary algorithm
The algorithm is based on the [evolution strategy](https://en.wikipedia.org/wiki/Evolution_strategy) framework, so there is only three rules:
1. change all current solutions by adding a normally distributed random vector (variation operator).
2. merge the changed solutions with the original ones (offspring recruitment).
3. select only the set of best solutions (elitism).

Solutions are encoded in the integer format (i.e., fixed searching grid). 
There is two options :
1. coarse grid: `uint8` 8 bit encoding, varying from 0 to 255.
2. fine grid: `uint16` 16 bit encoding, varying from 0 to 65535.

Solution decoding from gene to fenotype is performed by the following formula:
```markdown
f = (g * (f_hi - f_lo) / g_hi) + f_lo
```
Where `f` in the fenotype (real number), `g` is the genotype (integer number), 
`f_lo` and `f_hi` are the lower and upper bounds of the searching range of the fenotype (real numbers), respectively; 
and `g_hi` is the upper value of the searching grid (255 or 65535).

The strength of the variation operator is defined by the standard deviation `STD` 
of the normally distributed random vector. The vector is computed in the genetype domain (i.e., in integer values).
So the `STD` is computed from a rate parameter:
```markdown
STD = integer(R_STD * N_GRID)
```
Where `R_STD` is a rate value (0 to 1) and `N_GRID` is the grid resolution, so `N_GRID = g_hi`.

The algorithm is designed both for optimization and exploration.
When `EXPLORE = TRUE`, the procedure is set to find solutions within a certain range of the fitness score.
When `EXPLORE = FALSE`, the procedure is set to find solution that maximize (or minimize) the fitness score.

The pseudo-code:
```markdown
start
    get EXPLORE as boolean
    get N_POPSIZE as integer
    get N_GENERATIONS as integer
    get N_GRID as integer
    get R_STD as real
    set STD = R_STD * N_GRID as integer
    generate PARENTS as integer with POPSIZE
    evaluate PARENTS using FITNESS() function
    set g = 0
    repeat until g > GENERATIONS:
        set DELTA = NORMAL(mean=0, standard=STD) as integer
        set OFFSPRING = PARENTS + DELTA
        evaluate OFFSPRING using FITNESS() function
        merge OFFSPRING with PARENTS in POOL
        if EXPLORE == TRUE:
            sort POOL by exploration criteria
        else:
            sort POOL by fitness score
        select from POOL the best POPSIZE SAMPLE
        set PARENTS = SAMPLE
        set g = g + 1
    return PARENTS
end
```

Spoilers of the package - exploration in the Himmelblaus function: 

![anim](https://github.com/ipo-exe/evolvelab/blob/main/docs/spoiler.gif "spoiler")

### benchmark 2d functions

#### the `paraboloid` function

Equation:

```markdown
f(x, y) = level - (square(x - x0) + square(y - y0))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/parab.png "parab")

#### the `rastrigin 2d` function

Equation:

```markdown
f(x, y) = level - (20 + (square(x - x0) - 10 * cos(2 * pi * (x - x0))) + (square(y - y0) - 10 * cos(2 * pi * (y - y0))))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![rastr](https://github.com/ipo-exe/evolvelab/blob/main/docs/rastr_2d.png "rastr_2d")

#### the `himmelblaus` function

Equation:

```markdown
f(x, y) = level - (square(square(x - x0) + (y - y0) - 11) + square((x - x0) + square(y - y0) - 7))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=1000`, `x0=5`, `y0=5`:

![himm](https://github.com/ipo-exe/evolvelab/blob/main/docs/himm.png "himm")

### the `griewank` function

Equation:

```markdown
f(x, y) = level - 100 * (((square(x) + square(y)) / 4000) - (cos(x) * cos(y / sqrt(2))) + 1)
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/grie.png "grie")


## Results of fitting the 2d benchmark functions

Some experimental outputs for 2D (two decision variables). 

### Himmelblaus
Parameters: `N_GEN=100`, `N_POPSIZE=100`, `R_STD=0.7`

![himm_fit](https://github.com/ipo-exe/evolvelab/blob/main/docs/himm_fitting.gif "himm fitting")

Convergence plot:

![himm_conv](https://github.com/ipo-exe/evolvelab/blob/main/docs/himm_conv.png "himm convergence")

Retrieved scattergram (useful for uncertainty estimation)

![himm_scatter](https://github.com/ipo-exe/evolvelab/blob/main/docs/himm_scatter.png "himm scatter")

### Rastrigin
Parameters: `N_GEN=100`, `N_POPSIZE=200`, `R_STD=0.4`

![rastr_fit](https://github.com/ipo-exe/evolvelab/blob/main/docs/rastr_fitting.gif "rast fitting")

Convergence plot:

![rastr_conv](https://github.com/ipo-exe/evolvelab/blob/main/docs/rastr_conv.png "rast convergence")

Retrieved scattergram (useful for uncertainty estimation)

![rastr_scatter](https://github.com/ipo-exe/evolvelab/blob/main/docs/rastr_scatter.png "rast scatter")

### Parabola
Retrieved scattergram (useful for uncertainty estimation)

![parab_scatter](https://github.com/ipo-exe/evolvelab/blob/main/docs/parab_scatter.png "parab scatter")

### Griewank
Retrieved scattergram (useful for uncertainty estimation)

![grie_scatter](https://github.com/ipo-exe/evolvelab/blob/main/docs/grie_scatter.png "grie scatter")



## References

Beyer, HG., Schwefel, HP. Evolution strategies – A comprehensive introduction. Natural Computing 1, 3–52 (2002). https://doi.org/10.1023/A:1015059928466