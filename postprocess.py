import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

def create_outs(inputf, solf):
    """This function takes in as inputs the input file and the solution file
       and creates the pipe stencil, calculates mean temperature in the pipe,
       and plots the temperature distribution and overlays the mean temp. """

    #Reads files to get dimensions and spacing
    with open(inputf, 'r') as f:
        length, width, h = (float(x) for x in f.readline().split())

    #Read solution and find the length
    soln = np.loadtxt(solf)
    sol_len = len(soln)
    #Calculate x and y length
    xlen = int(length/h) + 1
    ylen = int(width/h) + 1

    #Create the pipe as an np 2D array 
    pipe = np.zeros((ylen, xlen), dtype = np.float64)
    #Add the hot isothermal boundary
    pipe[0] = soln[0:xlen]
    #Add the left periodic boundary 
    pipe[1:ylen-1, 0] = soln[xlen:sol_len-xlen:xlen]
    #Add the interior points and right periodic boundary
    pipe[1:ylen-1, 1:xlen] = soln[xlen:sol_len-xlen].reshape(ylen-2, xlen-1)
    #Add the cold isothermal boundary
    pipe[ylen-1] = soln[sol_len-xlen:sol_len]

    #Compute and output mean temperature
    meanT = np.mean(pipe)
    print("Mean Temperature: {:.5f}".format(meanT))

    #Create the plot
    plt.figure()
    plt.axis('equal')
    #Setup coordinate system
    x = np.arange(0, length + h, h)
    y = np.arange(0, width + h, h)
    X, Y = np.meshgrid(x, y)
    #Create colored plot with colorbar
    plt.pcolormesh(X,Y,pipe)
    #Add isoline and colorbar
    plt.colorbar()
    plt.contour(X, Y, pipe, [meanT])
    #Format plot
    plt.xlabel('x')
    plt.ylabel('y')
    #Save figure to see 
    plt.savefig("plot.png")

    return

if __name__ == "__main__":
    #Prompt usage
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python3 {} [input file] [soln file]". format(sys.argv[0]))
        sys.exit(0)

    #Initialize arguments
    inputf = sys.argv[1]
    solf = sys.argv[2]

    #Output statement
    print("Input file processed: {}".format(inputf))

    #Create remaining outputs
    create_outs(inputf, solf)

    
