# evolvelab
Evolutionary computing lab. For evolutionary computing stuff and experiments.

Spoilers of the package - evolution in the Himmelblaus function: 

![anim](https://github.com/ipo-exe/evolvelab/blob/main/docs/spoiler.gif "spoiler")


## 2d evolutionary algorithm

The pseudo-code:
```markdown
start
    get POPSIZE
    get GENERATIONS
    get ELITISM
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
level - (np.square(x - x0) + np.square(y - y0))
```
where `level`, `x0` and `y0` are pre-set parameters

Example of `level=100`, `x0=5`, `y0=5`:

![parab](https://github.com/ipo-exe/evolvelab/blob/main/docs/parab.png "parab")