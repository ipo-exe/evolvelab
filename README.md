# evolvelab
Evolutionary computing lab. For evolutionary computing stuff and experiments.

## The evolutionary algorithm
The algorithm is designed both for optimization and exploration.
When `EXPLORE = TRUE`, the procedure is set to find solutions within a certain range of the fitness score.
When `EXPLORE = FALSE`, the procedure is set to find solution that maximize (or minimize) the fitness score.

It is based on the [evolution strategy](https://en.wikipedia.org/wiki/Evolution_strategy) framework, so there is only three rules:
1.  change all solutions by adding a normally distributed random vector (variation operator).
2.  merge the changed solutions with the original ones (offspring recruitment).
3.  select only the set of best solutions (elitism).

The pseudo-code:
```markdown
start
    get EXPLORE as boolean
    get N_POPSIZE as integer
    get N_GENERATIONS as integer
    get N_GRID as integer
    get R_STD as real
    set STD = R_STD * N_GRID
    generate PARENTS with POPSIZE
    evaluate PARENTS using FITNESS() function
    set g = 0
    repeat until g > GENERATIONS:
        set DELTA = NORMAL(mean=0, standard=STD)
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
level - (square(x - x0) + square(y - y0))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/parab.png "parab")

#### the `rastrigin 2d` function

Equation:

```markdown
level - (20 + (square(x - x0) - 10 * cos(2 * pi * (x - x0))) + (square(y - y0) - 10 * cos(2 * pi * (y - y0))))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![rastr](https://github.com/ipo-exe/evolvelab/blob/main/docs/rastr_2d.png "rastr_2d")

#### the `himmelblaus` function

Equation:

```markdown
level - (square(square(x - x0) + (y - y0) - 11) + square((x - x0) + square(y - y0) - 7))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=1000`, `x0=5`, `y0=5`:

![himm](https://github.com/ipo-exe/evolvelab/blob/main/docs/himm.png "himm")

### the `griewank` function

Equation:

```markdown
level - 100 * (((square(x) + square(y)) / 4000) - (cos(x) * cos(y / sqrt(2))) + 1)
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/grie.png "grie")


##References



Beyer, HG., Schwefel, HP. Evolution strategies – A comprehensive introduction. Natural Computing 1, 3–52 (2002). https://doi.org/10.1023/A:1015059928466