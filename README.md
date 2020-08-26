# ENGR 491 - Computational Fluid Dynamics

Joshua Brinkerhoff

UBC Okanagan

joshua.brinkerhoff@ubc.ca

## How to use this resource

These exercises are inspired by Dr. Lorena Barba's [CFD Python: 12 Steps to Navier-Stokes](https://github.com/barbagroup/CFDPython/blob/master/README.md). They are intended to accompany Dr. Brinkerhoff's ENGR 491 CFD course at UBC Okanagan. The exercises roughly mirror the order that material is presented in ENGR 491. If students conduct about one exercise per week, they will be remain abreast of the lectures from ENGR 491.

Here is how I recommend students to use these exercises. 

1. Use them to see how the lecture concepts are implemented in a simple numerical solution. Sometimes the concept seems simple but the implementation in a real computer code can be tricky. The exercises provide simple examples that walk students through the implementation of core CFD concepts. Students will be able to see how the concepts can be implemented.

2. The philosophy behind these exercises is "learning by doing" and "learning by following". Students will learn by doing their work, but at the same time, complete examples are provided that students can follow. Combining both approaches will hopefully yield a very teachable and helpful way to learn. Therefore, students are encouraged to work through each exercise in an incremental way. Try to re-code each exercise on their own. Resist the urge to merely copy my solution--you will learn much more by wrestling through the challenges and bottlenecks that you're sure to encounter. At the same time, a complete solution of the working code is provided for you to follow.

3. Students are encouraged to experiment with their solutions or the solutions provided. Some suggestions for further learning are provided at the end of each exercise. Think through the questions that are presented there and try to use your code to explore your solution further. Explore, change parameters, and see how the solution behaves. You will gain an appreciation for the strengths and weaknesses of different methods, the sensitivity of a given solution, and begin to learn good programming practices through this kind of experimentation.

## Exercises

[Exercise 0 - Introduction to Python](https://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/00_Quick_Python_Intro.ipynb) This exercise is borrowed from Dr. Lorena Barba's excellent resource _12 Steps to Navier-Stokes_. Source: Barba, Lorena A., and Forsyth, Gilbert F. (2018). CFD Python: the 12 steps to Navier-Stokes equations. Journal of Open Source Education, 1(9), 21, [https://doi.org/10.21105/jose.00021](https://doi.org/10.21105/jose.00021)

[Exercise 1 - Numerical differentiation using finite differences](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/01_Exercise1.ipynb)

[Exercise 2 - Solving the steady 2D heat equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/02_Exercise2.ipynb)

[Exercise 3 - Solving the unsteady 2D heat equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/03_Exercise3.ipynb)

[Exercise 4 - Solving the 1D linear convection equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/04_Exercise4.ipynb)

[Exercise 5 - Solving the 1D non-linear convection equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/05_Exercise5.ipynb)

[Exercise 6 - Solving the 2D inviscid Burgers equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/06_Exercise6.ipynb)

[Exercise 7 - Solving the 1D Euler equation](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/07_Exercise7.ipynb)

[Exercise 8 - Solving the 2D Navier-Stokes equation in a cavity](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/08_Exercise8.ipynb)

[Exercise 9 - Fully-developed flow in a turbulent channel](https://nbviewer.jupyter.org/github/okcfdlab/engr491/blob/master/exercises/09_Exercise9.ipynb)

## Setup instructions

To use these lessons, you need Python 3, and the standard stack of scientific Python: NumPy, Matplotlib, SciPy, Sympy. And of course, you need [Jupyter](http://jupyter.org)—an interactive computational environment that runs on a web browser. The exercises are built as a set of [Jupyter notebooks](https://jupyter-notebook.readthedocs.org/en/latest/notebook.html) containing written introductions and explanations of the solution as well as complete working solutions implemented as Python code. To work through the exercises, start a new notebook and follow along with the provided materials, avoiding the temptation to copy and paste. As you face hurdles, the solution is available for you to consult. Once you have completed an exercise, test it to make sure you get the same output as what's provided. Experiment with your solution by changing solution parameters and see how they affect the accuracy and time required for the solution.

#### Installing via Anaconda
We *highly* recommend that you install the [Anaconda Python Distribution](http://docs.continuum.io/anaconda/install). It will make your life so much easier. 
You can download and install Anaconda on Windows, OSX and Linux. 

After installing, to ensure that your packages are up to date, run the following commands in a terminal:

```Bash
conda update conda
conda update jupyter numpy sympy scipy matplotlib
```

If you are unfamiliar with running terminal commands, you can use the Anaconda Navigator to ensure that the above packages are installed into your base environment.

#### Without Anaconda
If you already have Python installed on your machine, you can install Jupyter using pip:

```Bash
pip install jupyter
```

Please also make sure that you have the necessary libraries installed by running

```Bash
pip install numpy scipy sympy matplotlib
```

### Running the notebook server

Once Jupyter is installed, open up a terminal and then run 

```Bash
jupyter notebook
```

This will start up a Jupyter session in your browser!


