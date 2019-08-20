# Note: this does not conserve mass

import sys
import numpy as np
from scipy.stats import multivariate_normal

# SEED = 94550

def generate_gaussian(length, height, depth=None): 
    if depth == None: 
        depth = length
    x, y = np.meshgrid(np.linspace(-1, 1, length), np.linspace(-1, 1, depth))
    xy = np.column_stack([x.flat, y.flat])
    mu = np.array([0, 0])
    sigma = np.array([0.5, 0.5])
    covariance = np.diag(sigma**2)
    z = multivariate_normal.pdf(xy, mean=mu, cov=covariance)
    z = z.reshape(x.shape)
    return (height*z)

def generate(filename, Length, Height, Depth, gaussian_length, gaussian_height, gaussian_depth=None): 
    if gaussian_depth == None: 
        gaussian_depth = gaussian_length
        
    assert gaussian_height >= 2
    assert gaussian_height <= Height // 3
    assert gaussian_length >= gaussian_height
    assert gaussian_length <= 5 * gaussian_height
    assert gaussian_depth >= gaussian_height
    assert gaussian_depth <= 5 * gaussian_height
    
    num_gaussians = max(2 * 5 * Length * Depth // (gaussian_height * gaussian_length * gaussian_depth), 1)
    world = np.zeros((Length, Depth), dtype=int)
    for _ in range(num_gaussians): 
        pile = generate_gaussian(gaussian_length, gaussian_height, gaussian_depth)
        pile_start_i = np.random.randint(Length)
        pile_start_k = np.random.randint(Depth)
        for pile_i in range(gaussian_length): 
            for pile_k in range(gaussian_depth): 
                world_i = pile_start_i + pile_i
                if world_i >= Length: 
                    continue #out of bounds
                world_k = pile_start_k + pile_k
                if world_k >= Depth: 
                    continue #out of bounds
                world[world_i][world_k] += pile[pile_i][pile_k]
    world = world.astype(int)
    world = world.T
    np.savetxt(filename, world, fmt='%d')
    return world
                   
def main():
    #np.random.seed(SEED)
    if len(sys.argv) - 1== 6: 
        generate(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
    elif len(sys.argv) - 1 == 7: 
        generate(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]))
    else: 
        print("Need either 6 or 7 arguements: <filename>, <Length>, <Height>, <Depth>, <gaussian_length>, <gaussian_height>, <gaussian_depth>=gaussian_length \n")
    
        
main()

    
"""
from matplotlib.pyplot import imshow
world = generate("", 25, 25, 25, 7, 3)
print(world)
imshow(world)
"""
    
    
    
    
    
