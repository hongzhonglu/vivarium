# Composites

Cell models made of one or more processes that can be added to an environmental simulation.

# Execution

Add processs into an environmental simulation using ```experiment``` from agent framework. This is the same way the whole cell model is added to an environmental simulation: 

    python -m lens.environment.boot experiment --number N --type T
        
Here, ```T``` specifies the process type, such as ```chemotaxis```, and ```N``` specifies the number of cells.

Timelines are an optional argument:

    python -m lens.environment.boot experiment --number N --type T --timeline L

Here, ```L``` would be a string specifying events with time (seconds) and media_id: ```'0 minimal, 100 minimal_plus_amino_acids, 200 minimal'```

They can also be added to an already-running experiment with ```add```:

    python -m lens.environment.boot add --type T

    
# Making New Processes

# Making New Composites