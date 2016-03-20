/* Author: Phillip Levin
 * Date last modified: 03/19/2016
 * Program name: Prediction
 * Description: This is a C program that can generate quantitative predictions
 * based on a set of data. Users supply a set of x and y pairs, where x is a
 * suspected independent variable and y is a suspected dependent variable. After
 * naming the variables and providing a data set, users can provide a value and
 * the program will predict the corresponding dependent value. The algorithm
 * employed uses the supplied data set to determine a linear equation that
 * describes the relationship between the two variables. The program uses the
 * given pairs of values to generate a least squares equation by solving the
 * normal equation for the vector containing the regression coefficients that
 * define this model. Given this approach, Prediction is limited to making
 * predictions that are linear by nature. As a result, this program will not
 * produce optimal predictions for all data sets
 * Progress of development: Finished. Additional error checking and features me
 * be added in the future
 */

 /* Headers */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>

/* Constants */
#define PREDICT "Provide a value for %s, and I'll predict a value for %s. Or if
you are finished and want to end the program type !"
#define PREDICTION "I predict that %s will be %f"

/* normalEquationVariables structs are used to store information about the
 * variables b0 and b1, which are the coefficients in the equation y = b0 + b1x
 */
struct normalEquationVariables{

	double b0;
	double b1;
	int xCount;
	int yCount;
};

/***** Function definitions *****/

/* Function name: solveNormalEquation()
 * Function prototype:
 * void solveNormalEquation(struct normalEquationVariables
 *                          *v, double *xArr, double *yArr, int size)
 * Description: This function solves the normal equation XT*X*B=XTy, where X
 * is a matrix whose columns include the provided x values, XT is the inverse
 * of X, B is a vector containing the regression coefficients, and y is is a
 * vector comprising the y values. This function assumes that X is a 4x4 matrix
 * and that y is a 4x1 vector. It uses various linear algebra methods to solve
 * for the vector B and thus the regression coefficients. The
 * normalEquationVariables struct v is then updated with these values
 * Parameters:
 *				- v -- Struct whose members store the values of the regression
 *               coefficients b0 and b1
 *				- xArr -- Array containing supplied x values
 *				- yArr -- Array containing supplied y values
 *				- size  -- The number of elements in xArr and yArr
 *
 * Cause of error: No error checking performed
 * Return value: None
 */

void solveNormalEquation(struct normalEquationVariables *v, double *xArr,
                         double *yArr, int size)
{

	/* Make matrix X */
	double X[size][2];

	/* Each row of column 2 of X is determined by the values in xArr */
	int k; // row
	int j; // column
	for(k = 0; k < size; ++k)
	{
		X[k][0] = 1; // Coefficent of b0 is always 1
		X[k][1] = xArr[k]; // Coefficient of b1 is some x
	}

	/* Make matrix XT, which is the transpose of X */
	double XT[2][size];

	/* Each column of row 1 always has value of 1 */
	for(k = 0; k < size; ++k)
	{
		XT[0][k] = 1;
	}

	/* Each column of row 2 has the same value as the position in the same number
     * row of column 2 of matrix A
     */
	for(k = 0; k < size; ++k )
	{
		XT[1][k] = xArr[k];
	}

	/* Compute X^TX */

	/* XTX holds results of computation, which will be a matrix */
	double XTX[2][2];

	/* Compute column 1 */
	/* Compute colum 1, row 1 of X^TX */
	XTX[0][0] = size; // First row is always equal to the size of elements

	/* Compute column 1, row 2 of X^TX */
	double sum = 0; // Holds temporary sum of values during process matrix
                  // multiplication

	for(k = 0; k < size; ++k)
	{
		sum = sum + XT[1][k]; // Sum the values in each column of row 2
	}

	XTX[1][0] = sum;

	/* Compute column 2 */
	/* Compute column 2, row 1 of X^TX */
	XTX[0][1] = sum; // Matrix is symmetric

	/* Compute column 2, row 2 */
	sum = 0;

	for(k = 0; k < size; ++k)
	{
		sum = sum + X[k][1]*X[k][1]; // Sum the squares of each row of column 2
                                // of matrix X
	}

	XTX[1][1] = sum;

	/* Compute XTy */
	double XTy[2]; // Will always be two rows
	sum = 0;

	/* Compute row 1 of XTy */
	for(k = 0; k < size; ++k)
	{
		sum = sum + yArr[k]; // Sum each value in arrY
	}

	XTy[0] = sum;

	/* Compute row 2 of XTy */
	double xySum = 0; // Stores value of x*y for each pair of (x,y)
	sum = 0;

	for(k = 0; k < size; ++k)
	{
		xySum = xArr[k]*yArr[k]; // x*y for each pair of (x,y)
		sum = sum + xySum; // Sum each multiplication of pairs
	}

	XTy[1] = sum;

	/* Solve for B = ((XT*X)^-1)(y) using following steps */

	/* Solve for inverse of XT*X using */
	/* Need to find ad-bc in context of matrix XT*X */

	double ad = XTX[0][0] * XTX[1][1];
	double bc = XTX[0][1] * XTX[0][1];
	sum = ad-bc;
	sum = 1/sum; // 1/(ad-bc)

	/* Form matrix dbca /*/
	double dbca[2][2] = {
							XTX[1][1], -1*XTX[0][1], -1*XTX[1][0], XTX[0][0]

						};

	/* Multiply matrix adbc by 1/(ad-bc) */
	double inverseXTX[2][2];

	for(k = 0; k < 2; ++k)
	{
		for(j = 0; j < 2; ++j)
		{
			inverseXTX[k][j] = sum*dbca[k][j];
		}
	}

	/* Now to solve for B by multiplying the matrices ((XT*X)^-1)(y) */
	double vectorB[2];

	/* Solve for row 1 of Vector B */
	vectorB[0] = XTy[0]*inverseXTX[0][0] + XTy[1]*inverseXTX[0][1];
	vectorB[1] = XTy[0]*inverseXTX[1][0] + XTy[1]*inverseXTX[1][1];

	v->b0 = vectorB[0];
	v->b1 = vectorB[1];
}


/***** Main program *****/

/* Function name: main()
 * Function prototype: int main(int argc, char **argv)
 * Description: Main driver for program. Begins by prompting the user for the
 * name of the independent variable in question as well as a collection of data
 * points for that variable. Next, it prompts the user for the name of the
 * dependent variable and a set of data points that relate to each of the
 * previous x values provided. In the event that these data sets are not the
 * same size, it returns an error message to the user via stdout. If the sets
 * are of the same size, then it prompts user for an x value for which it will
 * predict a corresponding y value. It then calls solveNormalEquation() to
 * generate the coefficients of the linear equation describing the data. The y
 * variable in this equation is then solved for and printed to stdout. Users
 * are then prompted to enter another x value or end the program. Minimal error
 * checking occurs during execution of main() as well as the functions it calls.
 * Users are expected to follow the given directions, otherwise unexpected
 * program behavior may occur
 * Parameters:
 *            - Not used
 * Cause of error: If the size of the data sets for x and y variables are
 *                 different
 * Return value: 0 if program fails or ends early. 1 if program is run until
 *               completion
 */

int main(int argc, char **argv)
{
	/* Variables related to program functionality declared/defined */
	double x = 0; // Variable provided by user
	char xName[BUFSIZ]; // Name of x variable provided by user
	char yName[BUFSIZ]; // Name of y variable provided by user
	char *str = NULL; // Used as buffer for strings
	struct normalEquationVariables v; // Stores coefficients for variables b0
                                    // and b1 in linear equation y = b0 + b1x
	char temp[BUFSIZ]; // Line of text supplied by user is temporarily stored here
	double tempValue = 0; // Holds temporary double for each char provided by user
	double xArr[BUFSIZ]; // Holds each x value supplied by user
	double yArr[BUFSIZ]; // Holds each y value supplied by user

	/* Prompt user for data set for x values */
	printf("%s", "What is the name of your independent variable? Type ? if you
         need directions or ! at any time to end the program \n" );
	str = fgets(xName, BUFSIZ, stdin);
	if(*str == '!')
	{
		return 0;
	}

	if(*str == '?')
	{
		printf("%s", "Suppose you have two things that are related -- for example
    \nheight (independent variable) and shoe size (dependent variable). You
    \nthen collect data for a number of people, including their height and
    \nshoe size. Using this data, it is possible to then develop an equation
    \nthat will predict someone's shoe size based on their height. This
    \nprediction may or may not be completely accurate, depending on the
    true\nquantitative relationship between the two things, of course. Now to
    \nthe point of this program: using it, you can provide data for two
    \nrelated things and this the program will predict the value of one
    based\non the value of the other. So, let's try this again... \n");
		printf("%s", "What is the name of your independent variable? \n" );

		if(*fgets(xName, BUFSIZ, stdin) == '!'){
			return 0;
		}

	}

	/* Remove newline character at end of xName */
	str = strchr(xName, '\n');
	*str = '\0';

	printf("%s", "Please provide a set of data, with one entry per line, for ");
	printf("%s", xName );
	printf("%s", ". Enter a . on the last line to continue");
	printf("%c", '\n');

	int k = 0; // Keeps track of number of x values supplied by user

	while(*fgets(temp, BUFSIZ, stdin) != '.')
	{
		tempValue= atof(temp);
		xArr[k] = tempValue;

		++k;
	}

	v.xCount = k; // Update count for number of x values supplied

	/* Prompt user for data set for y value */
	printf("%s", "What is the name of your dependent variable? \n" );
	if(*fgets(yName, BUFSIZ, stdin) == '!')
	{
		return 0;
	}

	/* Remove newline character in yName */
	str = strchr(yName, '\n');
	*str = '\0';

	printf("%s", "Please provide a set of data, with one entry per line, for ");
	printf("%s", yName );
	printf("%s", ". Enter a . on the last line to continue");
	printf("%c", '\n');
	k = 0; // Keeps track of y values

	while(*fgets(temp, BUFSIZ, stdin) != '.')
	{
		tempValue= atof(temp);
		yArr[k] = tempValue;

		++k;
	}

	v.yCount = k; // Update count for number of y values supplied

	if(v.xCount != v.yCount)
	{
		printf("%s", "The number of x values must match the number of y values");
		return 0;
	}

	/* Prompt user for a value of x for which a y value will be predicted */
	printf(PREDICT, xName, yName);
	printf("%c", '\n');

	while(*fgets(temp, BUFSIZ, stdin) != '!')
	{

		x = atof(temp);

		/* Call solveNormalEquation to generate the values of b0 and b1 in the
         * linear equation of the form y = b0 + b1x
         */
		solveNormalEquation(&v, xArr, yArr, (int) k);

		/* Compute predicted value using linear equation */
		double y = 0;
		y = v.b0 + (v.b1*x);

		/* Print results */
		printf(PREDICTION, yName, y);
		printf("%c", '\n');

		printf(PREDICT, xName, yName);
		printf("%c", '\n');
	}

	return 1;
}
