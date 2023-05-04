# Exercise submissions for numerical relativity
### Prerequisites

Please install ```cmake``` by running

```
sudo apt-get update
sudo apt-get install cmake
```

### Usage
This project uses ```cmake``` to function, to build any particular executeable first navigate to the root directory of
the project and run

```
cmake ./
```

After the program has finished, build any target by running

```
make <target_name>
```

and run it by navigating to the output folder and executing the script

```
cd output/
./<target_name>
```

The target names and their functions for each project are found below.

### Targets
#### Exercise 1

* ```exercise1A``` - Tasks A.1, A.2 - Script for demonstrations of numerical approximations of the first and second derivatives of arbitrary functions
* ```exercise1B``` - Tasks B.1 through B.5 - Script for numerical computations involving the advection equation