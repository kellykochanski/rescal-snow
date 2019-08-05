__doc__ = '''utilities for creating, running, 
          saving, and visualizing for Conway\'s game of life'''

__author__ = 'Gian-Carlo DeFazio'


import numpy as np


a = np.array([[0,0,0,0,0],
              [0,1,1,1,0],
              [0,1,2,1,0],
              [0,1,1,1,0],
              [0,0,0,0,0]])

k = np.array([[0,0,0],
              [0,1,0],
              [0,0,0]])


k = np.array([[1,1,1],
              [1,0,1],
              [1,1,1]])
 


c = np.array([[0,0,0,0,0],
              [0,0,1,0,0],
              [0,0,1,0,0],
              [0,0,1,0,0],
              [0,0,0,0,0]])


g = np.array([[0,0,0,0,0,0,0,0],
              [0,0,1,0,0,0,0,0],
              [0,0,0,1,0,0,0,0],
              [0,1,1,1,0,0,0,0],
              [0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0]])


h = np.array([[0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0],
              [0,1,1,1,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,1,1,0,0],
              [0,0,0,0,0,0,1,1,0,0],
              [0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0]])

def convolve_2d(img, kernel):                                                                          
                                                                                                        
     # must deal with kernels with a dimension of 1 in either direction                                 
     # pad_vert                                                                                         
     p_v_l = kernel.shape[0] // 2                                                                       
     p_v_h = -(kernel.shape[0] // 2)                                                                    
     if kernel.shape[0] == 1:                                                                           
         p_v_h = img.shape[0]                                                                           
     #pad_horizontal                                                                                    
     p_h_l = kernel.shape[1] // 2                                                                       
     p_h_h = -(kernel.shape[1] // 2)                                                                    
     if kernel.shape[1] == 1:                                                                           
         p_h_h = img.shape[1]                                                                           
                                                                                                        
                                                                                                        
     v = img.shape[0]                                                                                   
     h = img.shape[1]                                                                                   
                                                                                                        
     img_filtered = None                                                                                
     if img.ndim == 3:                                                                                  
         img_filtered = np.zeros((img.shape[0] + kernel.shape[0] - 1,                                 
                                  img.shape[1] + kernel.shape[1] - 1,                                 
                                  img.shape[2]))                                                        
     else:                                                                                              
         img_filtered = np.zeros((img.shape[0] + kernel.shape[0] - 1,                                 
                                  img.shape[1] + kernel.shape[1] - 1))                                  
                                                                                                        
     # if factor of 0 do nothing                                                                        
     # if factor of 1 or -1, omit multiplication                                                        
     for (x,y), val in np.ndenumerate(kernel):                                                          
         if val == 1.0:                                                                                 
             img_filtered[x:x+v, y:y+h] += img                                                          
         elif val == -1.0:                                                                              
             img_filtered[x:x+v, y:y+h] -= img                                                          
         elif val:                                                                                      
             img_filtered[x:x+v, y:y+h] += (img * val)                                                  
                                                                                                        
     return img_filtered[p_v_l:p_v_h, p_h_l:p_h_h] 



def conv_con(img, kernel):
                                                                                                        
     # must deal with kernels with a dimension of 1 in either direction                                 
     # pad_vert                                                                                         
     p_v_l = kernel.shape[0] // 2                                                                       
     p_v_h = -(kernel.shape[0] // 2)                                                                    
     if kernel.shape[0] == 1:                                                                           
         p_v_h = img.shape[0]                                                                           
     #pad_horizontal                                                                                    
     p_h_l = kernel.shape[1] // 2                                                                       
     p_h_h = -(kernel.shape[1] // 2)                                                                    
     if kernel.shape[1] == 1:                                                                           
         p_h_h = img.shape[1]                                                                           
                                                                                                        
                                                                                                        
     v = img.shape[0]                                                                                   
     h = img.shape[1]                                                                                   
                                                                                                        
     img_filtered = None                                                                                

     img_filtered = np.zeros((img.shape[0] + kernel.shape[0] - 1,                                 
                              img.shape[1] + kernel.shape[1] - 1))                                  
                                                                                                        
     # if factor of 0 do nothing                                                                        
     # if factor of 1 or -1, omit multiplication                                                        
     for (x,y), val in np.ndenumerate(kernel):                                                          

         img_filtered[x:x+v, y:y+h] += img                                                          

                                                                                                        
     return img_filtered[p_v_l:p_v_h, p_h_l:p_h_h] 


 
# def convolve_2d(image, kernel):
#     '''expects kernel to be no lager than image in either dimension'''
#     # flip the kernel along both axes
#     kernel_flipped = np.flip(kernel, axis=(0,1))
    
    
#     # calculate size of output
#     h_image, w_image = image.shape
#     h_kernel, w_kernel = kernel.shape

#     h_output = h_image - h_kernel + 1
#     w_output = w_image - w_kernel + 1
#     output = np.zeros([h_image, w_image])

#     breakpoint()
#     for i in range(h_output):
#         for j in range(w_output):
#             output[i:i+h_kernel, j:j+w_kernel] += kernel_flipped * image[i:i+h_kernel, j:j+w_kernel]

#     return output[h_kernel//2:h_output+1, w_kernel//2:w_output+1]
    

def auto_convolve_2D(image):
    width = image.shape[0] // 2
    wrapped_image = np.pad(image, width, mode='wrap')
    breakpoint()
    return convolve_2d(wrapped_image, image)
    
    
# now make it stochastic
#def random_transitions(x)


def relu(x):
    maximum(x,0)



class Conway:
    '''This class implements the normal deterministic Conway's
    game of life. It also has options to make the transitions 
    probabilistic.'''

    
    def __init__(self, height=12, width=12, p=0.5, cells=None):
        if cells is None:
            temp_cells = np.random.random([height, width])
            # set the values to 1 and 0, set to 1 with probability p
            self.cells = np.empty([height, width])
            self.cells[temp_cells <= p] = 1
            self.cells[temp_cells > p] = 0
            self.cells = self.cells.astype(np.uint8)    
        else:
            self.cells = cells.astype(np.uint8)

        self.height, self.width = self.cells.shape
        
        self.kernel = np.array([[1,1,1],[1,0,1],[1,1,1]], dtype=np.uint8)

        
    def __repr__(self):
        return str(self.cells)


    def make_randos(self, p=0.5):
        allow_randos = np.random.random([self.height, self.width])
        return  (allow_randos <= p).astype(np.uint8)
        
    

    def transition(self):
        '''Do a normal, deterministic transition on the entire space.
        The edges wrap.'''
        wrapped = np.pad(self.cells, 1, mode='wrap')
        neighbor_counts = convolve_2d(wrapped, self.kernel)[1:-1,1:-1]
        self.cells = np.logical_or(neighbor_counts == 3,
                                   np.logical_and(neighbor_counts == 2, self.cells == 1)).astype(np.uint8)
        
        
        
    def transition_randomized(self, p_allow=0.99, p_spontaneous=0.002):
        '''Do a randomized transitions. For all cells that would change for a normal transition,
        transition with aprobability p_allow. Add a state change to each cell with probability 
        p_spontaneous. If a transition at a particular cell that would normally occur is prevented 
        from happening by p_allow, this will not be counter acted by p_spontaneous.
        That is, disallowing of transitions take precedence of spantaneous transitions at each cell.'''


        # if allow_randos is 1 at a particular cell, then normal transitions are allowed to occur
        allow_randos = np.random.random([self.height, self.width])
        allow_randos = (allow_randos <= p_allow).astype(np.uint8)

        # if spontaneous_randos is 1 at a particular cell, then a spontaneous transition will occur there
        spontaneous_randos = np.random.random([self.height, self.width])
        spontaneous_randos = (spontaneous_randos <= p_spontaneous).astype(np.uint8)
        spontaneous_randos = np.bitwise_and(allow_randos, spontaneous_randos)

        
        
        # do a normal transition, but don't overwrite current state
        wrapped = np.pad(self.cells, 1, mode='wrap')
        neighbor_counts = convolve_2d(wrapped, self.kernel)[1:-1,1:-1]
        next_states= np.logical_or(neighbor_counts == 3,
                                   np.logical_and(neighbor_counts == 2,self.cells == 1)).astype(np.uint8)

        current_states = self.cells.copy()
        
        # now only allow chnges where randos == 1 to happen
        # determine which cells changes will happen for

        
        self.cells = np.bitwise_xor(np.bitwise_and(np.bitwise_xor(self.cells, next_states), allow_randos), self.cells)

        # now add the spontaneous transitions
        self.cells = np.bitwise_xor(self.cells, spontaneous_randos)
        

    def transition_spontaneous(self, p=0.0015):
        randos = self.make_randos(p)
        self.transition()
        self.cells = np.bitwise_xor(self.cells, randos)




        
    def make_sequence(self, filename, n=20):
        '''create a sequence of n transitions and saves as a 3D array'''
        sequence = np.empty([n, self.height, self.width], dtype=np.uint8)
        sequence[0] = self.cells[:,:]
        
        for i in range(1,n):
            self.transition()
            sequence[i] = self.cells[:,:]
            print(sequence[i])
        np.savez_compressed(filename, sequence, 's')
            

        

if __name__ == '__main__':
    c = Conway(cells=g)
    print(c)
    for i in range(50):
        print('\n' + str(i) + '\n')
        c.transition_spontaneous()
        print(c)

