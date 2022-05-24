# evolvelab
Evolutionary computing lab. For evolutionary computing stuff and experiments.

Spoilers of the package - evolution in the Himmelblaus function: 

![anim](https://github.com/ipo-exe/evolvelab/blob/main/docs/spoiler.gif "spoiler")


## 2d evolutionary algorithm

The pseudo-code:
```markdown
start
    get N_POPSIZE
    get N_GENERATIONS
    generate PARENTS with POPSIZE
    evaluate PARENTS fitness scores
    set g = 0
    repeat until g > GENERATIONS:
        generate OFFSPRING using PARENTS
        evaluate OFFSPRING fitness scores
        merge OFFSPRING with PARENTS in POOL
        select next generation PARENTS from POOL
        if ELITISM == True:
            sort POOL by fitness scores
            select from POOL the best POPSIZE SAMPLE
        else:
            select from POOL a POPSIZE SAMPLE using RWS method
        set PARENTS = SAMPLE
        set g = g + 1
    return PARENTS
end
```

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

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/parab.png "grie")