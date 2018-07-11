import os
import sys
import time as t
import numpy as np
import pandas as pd

#Reads a file (1 row per line, each row should have same number of column values, separated by whitespaces), returns a numpy array
#filename -> the name of the file to open
#datatype -> the datatype of the values within the file (int, or float)
def read_data(filename,datatype):
    
    sig_2d = []
    f = open(filename, 'r')
    try:
        for line in f.readlines():
            line=line.strip()
            sig_2d.append(map(datatype,line.split()))
    finally:
        f.close()
        return np.array(sig_2d)

#Reads all files in directory with specific extension and creates a list of numpy arrays containing numerical data
#dir_path -> The directory containing data files
#ext -> The extension of datafiles
#datatype -> The type of data stored in each files (int or float data)
def read_directory(dir_path,ext,datatype):
    
    #Get array list of all files to open
    files = []
    try:
        for file in os.listdir(dir_path):
            if file.endswith(ext):
                files.append(os.path.join(dir_path,file))
    except:
        print("An error occured, check correct directory was passed.")    

    files.sort()
    count = 0.0
    max = len(files)
    #Create list of numpy arrays containing data for each file
    all_data = []
    for file in files:
        all_data.append(read_data(file,datatype))
        count += 1
        progress = (count / max) * 100.0
        sys.stdout.write('\r[{}] {}%'.format('#'*int(progress/5), round(progress,2)))
        sys.stdout.flush()
    return all_data


#Performs fft2d analysis on the data taken from input file and returns single quadrant fft result.
#data -> data to perform fft on
def fft2d_analyze(data):

    #Data points for x and y axis
    dpx, dpy = data.shape

    #Create x, y axis for graphing 2d
    x = np.arange(0,dpx)
    y = np.arange(0,dpy)

    #Set DC frequency bin to 0
    fft_data = data - np.mean(data)

    #Get fft2d and resize to single quadrant
    return np.fft.fft2(fft_data)[0:dpx/2,0:dpy/2]*2/(dpx*dpy)

#Performs fft2d analysis on all data in data list
#all_data -> a list of 2d numpy arrays containing data
def all_fft2d_analysis(all_data):
    
    all_fft_results = []
    for data in all_data:
        all_fft_results.append(fft2d_analyze(data))
    
    return all_fft_results

#Calculates all the amplitudes of the fft_analysis results
#all_fft -> fft results taken from all the data in a directory
def all_amplitudes(all_fft):
    
    all_amps = []
    for data in all_fft:
        all_amps.append(np.abs(data))
    return all_amps

#Locates the dominant frequencies within amplitude data by filtering with a threshold value
#Returns a list of dominant frequencies as well as a full list of all data of interest
#all_amps -> a list of numpy arrays containing amplitude data from fft results
#threshold -> the value used to filter out and detect dominant frequencies
def get_dominant_freqs(all_amps, threshold):
    dominant_freqs = []
    for amps in all_amps:
        top = top_percent(amps, threshold)
        if len(top) == 1:
            coords = (top[0][1],top[0][2])
            if coords not in dominant_freqs:
                dominant_freqs.append(coords)
    return dominant_freqs

#Get the top values based on a percentile threshold.
#E.g. .8 threshold will return data points with values >= (max-min)*.8 + min.
#values -> 2d array containing the values to threshold
#threshold -> thresholding percentage as a decimal (0.8=80%)
def top_percent(values, threshold):
    
    min = np.amin(values)
    max = np.amax(values)
    thresh_value = (max-min)*threshold+min
    
    top = []
    
    for i, y in enumerate(values):
        for ii, x in enumerate(y):
            if x >= thresh_value:
                top.append([x,ii,i])
                  
    return top

#Returns 
def top_values(values, threshold, d_freqs):
    
    min = np.amin(values)
    max = np.amax(values)
    thresh_value = (max-min)*threshold+min
    
    top = []
    
    for i, y in enumerate(values):
        for ii, x in enumerate(y):
            if (ii, i) in d_freqs:
                top.append([x,ii,i,True])
            elif x >= thresh_value:
                top.append([x,ii,i,False])
                  
    return top

#Returns the fft results as a Pandas dataframe 
#timestep -> the current frame/time of the data being passed
#threshold -> the threshold value of data to track when providing the stats
#amps -> 2d array that contains the amplitude results of fft
#phases -> 2d array that contains the phase results of fft
#d_freqs -> a list of the dominant frequencies found among all data
def get_stats(timestep,threshold,fft_data,amps,d_freqs):

    #get stats for top values
    top_vals = top_values(amps, threshold, d_freqs)

    stats = []

    for i, data in enumerate(top_vals):
        amp = data[0]
        x = data[1]
        y = data[2]
        dom = data[3]
        phase = np.angle(fft_data[y][x])
        if x > 0:
            wave = 1.0/x
        else:
            wave = 1
        aspect = amp/wave
        stats.append([timestep,dom,x,y,amp,phase,wave,aspect])
    
    return pd.DataFrame(stats,columns=['Time','Dominant Freq.','X','Y','Amplitude','Phase','Wavelength','Amp/Wave'])

#Returns one large dataframe containing all the results taken from fft analysis of all files in directory
#threshold -> the threshold value of data to track when providing the stats
#d_freqs -> a list of the dominant frequencies found among all data
#all_fft_data -> a list of numpy arrays containing the fft data taken from all the files
#all_amps -> a list of numpy arrays conatining all the amplitude data
def get_all_stats(threshold, d_freqs, all_fft_data, all_amps):
    
    frames = []
    for i, amp in enumerate(all_amps):
        frames.append(get_stats(i,threshold,all_fft_data[i],amp,d_freqs))

    #Create one large dataframe
    master_frame = pd.concat(frames[1:])
    
    return master_frame
    
#Graphs desired data at specific intervals and saves the graph as a png to specified directory, creates a GIF of all images and saves with specified name
#Returns the total number of images produced
#int -> the plot interval, e.g. 10 means every 10th element in the data list will be graphed and saved as png, 0 means no graphing
#dir -> the directory to save the graphs to
#GIF_name -> the full file name (directory + name + .gif) to use when creating the GIF animation
#basename -> the base filename to use, each name will be given a index number to make it unique, e.g "filename{}" -> filename0.png, filename1.png... etc.
#data -> a list of several numpy arrays of 2d data for plotting
#graph_type -> plot type to use, 'surf'->surface, 'wire'->wireframe, 'scat'->scatter, 'cont'->contour
#fig_size -> (optional) size of the figure to plot
#x_label,y_label,z_label -> (optional) labels to use for x,y and z-axis
#title -> (optional) the title to use for plot, note if you have a title like: 'title a {id}' the {id} will be replaced by the index of that data snapshot
def graph_all(interval,dir,GIF_name,basename,data,graph_type,fig_size,x_label,y_label,z_label,title):

    if interval == 0:
        return 0

    count = 0.0
    max = len(data)/interval
    images = []
    for i, d in enumerate(data):
        if i % interval == 0:
            count += 1
            progress = count / max * 100
            sys.stdout.write('\r[{}] {}%'.format('#'*int(progress/5), round(progress,2)))
            sys.stdout.flush()
            fig = plot_data(d,graph_type,fig_size,x_label,y_label,z_label,title.format(id=i))
            fname = dir + basename.format(i)
            save_to_png(fig,fname)
            pl.close()
            images.append(imageio.imread(fname+'.png'))
    
    #Create gif animation from images
    imageio.mimsave(GIF_name,images)

    return len(images)

#Creates a plot of amplitude data and returns the figure
#data -> a numpy array of 2d data for plotting             
#type -> plot type to use, 'surf'->surface, 'wire'->wireframe, 'scat'->scatter, 'cont'->contour
#fig_size -> (optional) size of the figure to plot
#x_label,y_label,z_label -> (optional) labels to use for x,y and z-axis
#title -> (optional) the title to use for plot
def plot_data(data,type='wire',fig_size=(10,10),x_label='x',y_label='y',z_label='Amplitude',title='Graph of FFT2d Amplitude Spectrum'):
    w, l = data.shape
    data_x, data_y = np.meshgrid(np.arange(0,w),np.arange(0,l), indexing='ij')        
    fig = pl.figure(figsize=fig_size)
    ax = fig.add_subplot(111, projection='3d')
    ax.set(xlabel=x_label,ylabel=y_label,zlabel=z_label,title=title)
    if type == 'wire':
        ax.plot_wireframe(data_x,data_y,data)
    elif type == 'surf':
        ax.plot_surface(data_x,data_y,data)
    elif type == 'scat':
        ax.scatter(data_x,data_y,data,s=5)
    elif type == 'cont':
        ax.contour(data_x,data_y,data)
    else:
        ax.plot_wireframe(data_x,data_y,data)

    return fig
    
#Saves the figure as a png image
def save_to_png(figure, fname):
    figure.savefig(fname+'.png',bbox_inches='tight')
    
#Performs fft2d analysis on the data taken from input file.
#time_step -> The time value of current data
#data -> a 2d numpy array containing numerical values
#input_type -> The data type of the input (int or float values for example)
#graph -> if true a figure with graphs will be made using pyplot and function will return the figure, otherwise 0
#both -> if true, both the amplitude and input data will be graphed (if graph is true), otherwise only amplitude data is graphed
def fft2d_analysis(time_step, threshold, data, graph, both):

    #Data points for x and y axis
    dpx, dpy = data.shape

    #Create x, y axis for graphing 2d
    x = np.arange(0,dpx)
    y = np.arange(0,dpy)

    #Create x-y plane for graphing 3d
    X, Y = np.meshgrid(x, y, indexing='ij')

    #Set DC frequency bin to 0
    file_data = data - np.mean(data)

    #Get fft2d and resize to single quadrant
    fft2 = np.fft.fft2(file_data)[0:dpx/2,0:dpy/2]*2/(dpx*dpy)
    
    #Get amplitude spectrum
    amp_s = np.abs(fft2)
    amp_x, amp_y = np.meshgrid(np.arange(0,dpx/2),np.arange(0,dpy/2), indexing='ij')    

    #Get phase 
    phases = np.angle(fft2)

    #Get the stats from data
    stats = get_stats(time_step,threshold,amp_s,phases)

    if graph:
        fig = plot_data([X, Y, data],[amp_x,amp_y,amp_s],both)

        return fig, stats

    return stats

#directory -> The directory that holds the log files to analyze
#image_interval -> A graph will be made and saved at every interval. E.g 50 = every 50th data file will be graphed.
#Note: if image interval is set to 0, no graphs are made.
def main(directory="input_data/ALT_DATA1/",output_dir="results_1",image_interval=100,base_ext='.data'):

    MAIN_DATA_DIR = directory
    PNG_OUTPUT_DIR = output_dir + "png_output/"
    DATA_OUTPUT_DIR = output_dir
    CSV_OUTPUT_NAME = "fft_analysis.csv"
    GIF_OUTPUT_NAME = "fft_results.gif"
    GRAPH_TYPE = 'surf' #Options available to use, 'surf'->surface, 'wire'->wireframe, 'scat'->scatter, 'cont'->contour
    XLABEL = 'x'
    YLABEL = 'y'
    ZLABEL = 'Amplitude'
    TITLE = 'Graph of FFT2d Amplitudes at t={id}' #id is replaced by index of snapshot
    FIG_SIZE = (10,10)
    SNAPSHOT_INTERVAL = image_interval #A png image of graph is saved at each interval
    BASE_FILE_NAME = "ALTI{:05d}_t0"
    BASE_EXT = base_ext
    THRESHOLD = 0.8
    DATA_TYPE = int

    #Check parent directories exists
    directories = [MAIN_DATA_DIR, DATA_OUTPUT_DIR]
    for d in directories:
        if not os.path.exists(d):
            print("The specified path: {} was not found. Analysis cancelled.".format(d))
            return 1

    #Make png output directory if needed
    if image_interval > 0 and not os.path.exists(PNG_OUTPUT_DIR):
        os.makedirs(PNG_OUTPUT_DIR)

    #Get all data from the main directory
    print("Reading data files...")
    t0 = t.time()
    all_data = read_directory(MAIN_DATA_DIR,BASE_EXT,DATA_TYPE)
    t1 = t.time()

    #Check input files exists
    file_count = len(all_data)
    if file_count == 0:
        print("No files with the extension: '{}' were found. Analysis cancelled.".format(BASE_EXT))
        return 1
    t_read = t1-t0
    print("\rTotal files read: {} (with ext: '{}')\nRead time: {}s\nCalculating fft2d data on all files...".format(file_count,BASE_EXT,t_read))
    
    #Get all fft2d data for calculating
    t0 = t.time()
    all_fft2d = all_fft2d_analysis(all_data)
    t1 = t.time()
    t_fft2d = t1-t0
    print("FFT calculations time: {}s\nCalculating amplitude data...".format(t_fft2d))

    #Get all amplitude data
    t0 = t.time()
    all_amps = all_amplitudes(all_fft2d)
    t1 = t.time()
    t_amps = t1-t0
    print("Amplitude calculation time: {}s\nUsing amplitude data to find dominant frequencies...".format(t_amps))

    #Find dominant frequencies
    t0 = t.time()
    d_freqs = get_dominant_freqs(all_amps,THRESHOLD)
    t1 = t.time()
    t_freqs = t1-t0
    print("{} dominant frequencies found at threshold {}%, time: {}s\nConcatenating all data and writing to csv...".format(len(d_freqs),THRESHOLD*100,t_freqs))
    
    #Concatenate and save data results
    t0 = t.time()
    all_stats = get_all_stats(THRESHOLD,d_freqs,all_fft2d,all_amps)
    all_stats.to_csv(DATA_OUTPUT_DIR + CSV_OUTPUT_NAME)
    t1 = t.time()
    t_stats = t1 - t0
    print("Analysis results complete time: {}s\nPlotting data at intervals of {} and creating PNG images...".format(t_stats,SNAPSHOT_INTERVAL))

    #Graph resulting data at specific intervals and save as images to directory, create GIF of pngs
    if image_interval > 0:
        #Import graphing libraries if needed
        import imageio
        import matplotlib.pyplot as pl
        from mpl_toolkits.mplot3d import Axes3D as pl3d

        t0 = t.time()
        im_count = graph_all(SNAPSHOT_INTERVAL,PNG_OUTPUT_DIR,DATA_OUTPUT_DIR+GIF_OUTPUT_NAME,BASE_FILE_NAME,all_amps,'surf',FIG_SIZE,XLABEL,YLABEL,ZLABEL,TITLE)
        t1 = t.time()
        t_plot = t1-t0
        t_total = t_read + t_fft2d + t_amps + t_freqs + t_stats + t_plot
        print("\r{} PNG's created in: {}s.\nGIF animation complete.\nAnalysis process complete!\nTotal time: {}s".format(im_count,t_plot,t_total))
    else:
        t_total = t_read + t_fft2d + t_amps + t_freqs + t_stats
        print("No PNG images or GIF animation made.\nAnalysis process complete!\nTotal time: {}s".format(t_total))
args = sys.argv

if len(args) > 4:
    main(args[1].args[2],int(args[3]),args[4])
elif len(args) > 3:
    main(args[1],args[2],int(args[3]))
elif len(args) > 2:
    main(args[1],args[2])
elif len(args) > 1:
    main(args[1])
else:
    main()
