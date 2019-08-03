import heightmap
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


# Data for plotting

def line_graph(datums):
    d = datums
    if not isinstance(datums, list):
        d = [datums]
    t = list(range(len(d[-1])))

    
    fig, ax = plt.subplots()

    for datum in d:
        ax.plot(t, datum)

    ax.set(xlabel='time (t_0)', ylabel='variance',
           title='who cares')
    ax.grid()

    plt.show()


def line_mods(datums):
    # expect 2d array

    t = list(range(datums.shape[1]))
    
    fig, ax = plt.subplots()

    for datum in datums:
        ax.plot(t, datum)

    ax.set(xlabel='time (t_0)', ylabel='variance',
           title='who cares')
    ax.grid()

    plt.show()


def line_all(md1, md2, title, filename, cd=None):
    # expect 2d array


    t = list(range(len(md1[0][0])))
    #t = list(range(len(cd)))
    
    fig, ax = plt.subplots(figsize=(24.5,16))
    ax.set_xlim(0.0, 17.)
    #ax.set_ylim(0.0, 3.0)
    ax.set_yscale('log')


    if cd is not None:
        ax.plot(t, cd, 'deeppink')
    for ri in md1:
        for mn in ri:
            ax.plot(t, mn, 'aqua', alpha=0.4)
    for ri in md2:
        for mn in ri:
            ax.plot(t, mn, 'blueviolet', alpha=0.25)


    # set the averages
    md1_ave = np.mean(md1, axis=(0,1))
    md2_ave = np.mean(md2, axis=(0,1))
    ax.plot(t, md1_ave, 'deepskyblue', alpha=1, linewidth=14)
#    ax.plot(t, md1_ave, 'white', alpha=1, linewidth=8)        
    ax.plot(t, md1_ave, 'aqua', alpha=1, linewidth=8, linestyle='-')
    ax.plot(t, md2_ave, 'indigo', alpha=1, linewidth=14)
#    ax.plot(t, md2_ave, 'white', alpha=1, linewidth=8)    
    ax.plot(t, md2_ave, 'blueviolet', alpha=1, linewidth=8, linestyle='-')

            
    legend_names = []
    if cd is not None:
        legend_names = ['sinusoidal', 'rectangular prism', 'control group']
    else:
        legend_names = ['sinusoidal', 'rectangular prism']
            
    ax.legend(legend_names, fontsize=30)
    leg = ax.get_legend()
    leg.legendHandles[0].set_color('aqua')
    leg.legendHandles[1].set_color('blueviolet')
    if cd is not None:
        leg.legendHandles[2].set_color('deeppink')


    
    ax.set_xlabel('time steps x 100', fontsize=26)
    ax.set_ylabel('average variance', fontsize=26)
    ax.tick_params(axis='both', which='major', labelsize=26)
    ax.tick_params(axis='both', which='minor', labelsize=28)    
    ax.set_title(title, fontsize=28)

    #ax.grid()
    plt.savefig(filename, transparent=True, bbox_inches='tight')

    plt.show()


    
# deal with control first
#


def c_mean(c):
    return np.mean(c, axis=(0,1))

# swap axes to iterate over rando
def c_mean_singles(c):
    cw = np.swapaxes(c,0,1)
    return np.mean(cw, axis=(0))

def c_variances_singles(c):
    cms = c_mean_singles(c)
    cw = np.swapaxes(c,0,1)
    return np.mean(np.square(cms - cw), axis=(0))

def control_variances(c):
    cm = c_mean(c)
    return np.mean(np.square(c - cm), axis=(0,1))

# for single experiment
def mod_mean_single(m):
    mw = np.swapaxes(m, 0,3)
    return np.mean(m,axis=0)

    
def mod_variances_singles(m):
    mms = mod_mean_single(m)
    mw = np.swapaxes(m, 0,3)
    return np.mean(np.square(mms - mw), axis=(0))




def make_pics(c, m, left=True, bottom_left=True):
    # left side
    if left:
        h0 = heightmap.HeightMap(c[4][4][0])
        h10 = heightmap.HeightMap(c[4][4][10])
        h30 = heightmap.HeightMap(c[4][4][30])
        h30.save_as_pdf('pics/h_4_4_30.pdf')
        plt.clf()
        h10.save_as_pdf('pics/h_4_4_10.pdf')
        plt.clf()        
        h0.save_as_pdf('pics/h_4_4_0.pdf')
        plt.clf()        

        h30.save_fft_blur('pics/f_4_4_30.pdf', in_3d=True)
        plt.clf()        
        h10.save_fft_blur('pics/f_4_4_10.pdf', in_3d=True)
        plt.clf()        
        h0.save_fft_blur('pics/f_4_4_0.pdf', in_3d=True)
        plt.clf()        

    # bottom left
    if bottom_left:
        b1 = heightmap.HeightMap(c[6][7][15])
        b2 = heightmap.HeightMap(c[3][2][15])
        b1.save_as_pdf('pics/h_6_7_15.pdf')
        plt.clf()        
        b2.save_as_pdf('pics/h_3_2_15.pdf')
        plt.clf()        
        b1.save_fft_blur('pics/f_6_7_15.pdf', in_3d=True)
        plt.clf()        
        b2.save_fft_blur('pics/f_3_2_15.pdf', in_3d=True)
        plt.clf()        
        b3 = heightmap.HeightMap(m[3,1,4,2,0])
        b4 = heightmap.HeightMap(m[5,0,5,6,0])
        b3.save_as_pdf('pics/h_3_1_4_2_4.pdf')
        plt.clf()        
        b4.save_as_pdf('pics/h_5_0_5_6_4.pdf')
        plt.clf()        
        b3.save_fft_blur('pics/f_3_1_4_2_4.pdf', in_3d=True)
        plt.clf()        
        b4.save_fft_blur('pics/f_5_0_5_6_4.pdf', in_3d=True)
        plt.clf()        

        
        # bottom middle
        # show progressions for b1
        for i in range(2,18,3):
            con = heightmap.HeightMap(c[6][7][i])
            con.save_as_pdf('pics/h_cont_' + str(i) + '.pdf')
            con.save_fft_blur('pics/f_cont_' + str(i) + '.pdf', in_3d=True)
            plt.clf()

            m0 = heightmap.HeightMap(m[3,1,4,2,i])
            m0.save_as_pdf('pics/h_mod0_' + str(i) + '.pdf')
            m0.save_fft_blur('pics/f_mod0_' + str(i) + '.pdf', in_3d=True)
            plt.clf()
            
            m1 = heightmap.HeightMap(m[5,0,5,6,i])
            m1.save_as_pdf('pics/h_mod1_' + str(i) + '.pdf')
            m1.save_fft_blur('pics/f_mod1_' + str(i) + '.pdf', in_3d=True)
            plt.clf()
        

def make_graphs(cf, mf):
    cm = c_mean(c)[:18]
    cv = control_variances(c)[:18]

    cv1 = np.mean(cv, axis=(1,2))

    # make a graph that looks at the variances within the subgroups
    get_all_datums(cf, mf)
    mv = mod_variances_singles(mf)
    mvv = np.mean(mv, axis=(4,5))
    line_all(mvv[0], mvv[1], 'Variances Within Modified Simulation Groups', 'pics/small_graph.pdf')
    
        
# get the control mean
# get the control variance
# get the mod variances from the control mean for each of 16 experiments
def get_all_datums(c, m):
    cm = c_mean(c)[:18]
    cv = control_variances(c)[:18]


    #mw = np.swapaxes(m , 0, 3)
    #print(mw.shape)
    v = np.square(m-cm)
    
    va = np.mean(v, axis=3)



    va1 = np.mean(va, axis=(4,5))
    cv1 = np.mean(cv, axis=(1,2))
    va1 = np.swapaxes(va1, 0, 1)
    print(va1[0].shape)
    print(cv1.shape)

    #va1 = np.mean(va1, axis=(2))
    
    line_all(va1[0],va1[1], 'Average Squared Difference of Modified Simulations Compared to Control Group Average',
             'pics/big_graph.pdf', cv1)
#    line_all(va1[1], cv1)
def hm(x):
    plt.imshow(x)
    plt.colorbar()
    plt.show()




if __name__ == '__main__':
    pass
    # get the files
    mf = np.load('mod_clipped_fft.npz')['arr_0']
    cf = np.load('control_clipped_fft.npz')['arr_0']
    m = np.load('mod_clipped.npz')['arr_0']
    c = np.load('control_clipped.npz')['arr_0']

    # #cv = np.mean(control_variances(cf), axis=(1,2))
    
    # #line_graph(cv)


    #make_pics(c, m)

    make_graphs(cf, mf)
    
    
    #cvs = np.mean(c_variances_singles(cf), axis=(2,3))
    #mvs = np.mean(mod_variances_singles(mf), axis=(4,5))





    
    #line_graph(cvs)
