Author: Phillip Levin

Name of program: Prediction:

Languages: C

-----------
Description
-----------
This is a C program that can generate quantitative predictions based on a set of data. Users supply a set of x and y pairs, where x is a suspected independent variable and y is a suspected dependent variable. After naming the variables and providing a data set, users can provide a value and the program will predict the corresponding dependent value. The algorithm employed uses the supplied data set to determine a linear equation that describes the relationship between the two variables. The program uses the given pairs of values to generate a least squares equation by solving the normal equation for the vector containing the regression coefficients that define this model. Given this approach, Prediction is limited to making predictions that are linear by nature. As a result, this program will not produce optimal predictions for all data sets.

---------------------
Status of development
---------------------
As of 03/19/2016, I am finished developing this program. The purpose of the program was to test out an idea I had for utilizing least squares/regression to develop models that could predict information for the user. The present program fufills this task to the degree I had envisioned for now. In the future, however, I may update it with some additional features as well as improved error checking. 

----------------------
How to compile program 
----------------------
Can be compiled using gcc with the command:

gcc main.c

------------------
How to run program
------------------
After compiling, the program can be run with the command:

./a.out
