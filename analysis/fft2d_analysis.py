"""
Rescal-snow: a cellular automaton model of self-organized snow
Copyright (C) 2019 Kelly Kochanski

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.
This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
import sys
sys.path.insert(0, "../scripts/utilities") #Find rescal_utilities here
import time as t
import numpy as np
import pandas as pd
import rescal_utilities as ru
import multiprocessing as mp
from scipy.optimize import curve_fit

can_plot = True
try:
    #Import graphing libraries if needed
    import imageio
    import matplotlib.pyplot as pl
    from mpl_toolkits.mplot3d import Axes3D as pl3d
    
except:
    can_plot = False
    pass

#Reads ReSCAL simulation output file (the input files for postprocessing)
#Reads a file (1 row per line, each row should have same number of column values, separated by whitespaces), returns a numpy array
#filename -> the name of the file to open
#datatype -> the datatype of the values within the file (int, or float)
def read_data(filename,datatype,show_error=True):
    sig_2d = []
    f = open(filename, 'r')
    try:
        for line in f.readlines():
            line=line.strip()
            # PY3 Note: map is iterable, numpy requires a list
            sig_2d.append([datatype(x) for x in line.split()])
            # sig_2d.append(map(datatype,line.split()))
    finally:
        f.close()
        arr = np.array(sig_2d)
        try:
            x, y = arr.shape
            # PY3 Note: indexing 0-50 because input data .log files are repetitive (imshow to verify)
            # return np.array(sig_2d)[0:50]
            return np.array(sig_2d)
        except:
            if show_error:
                print("Error in data in file: {}\nData is not formatted correctly or missing. Data was skipped.".format(filename))
            return None
        return np.array([])

#Writes to a textfile the basic summary of the data, using parameter file and summary dataframe
#run_stats = [input_directory, total_files, analysis_time, d_freq_count,dominant_threshold_percentage,amplitude_threshold]
def write_summary(run_stats,directory,filename,par_file_path,summary_data,top_amp_data):
    params = ru.Parameters()
    try:
        params.read(par_file_path)
    except:
        print("Parameter file not found at: {}".format(par_file_path))
        params = None

    filepath = directory+filename
    f = open(filepath,"w")

    try:
        #Function for wave_time curve_fit
        def func(t,a,b,c,d):
            return a*np.exp(-t/b)-c*t+d
        
        #Write summary data
        if params:
            in_dir = run_stats[0]
            f_count = run_stats[1]
            a_time = run_stats[2]
            d_freq_count = run_stats[3]
            d_thresh = run_stats[4]*100
            thresh = run_stats[5]
            s_time = top_amp_data["Time"].iloc[-1]
            
            freqs, fields = summary_data.shape
            lam_s = params.get("Lambda_S")
            tau = params.get("Tau_min")
            Ava = params.get("Ava_angle")
            h = params.get("H")
            d = params.get("D")
            l = params.get("L")
            f.write("Analysis Info:\nInput Directory: {}\nTotal files read: {}\nAnalysis time: {}\nSimulation Time: {}\nDominant Threshold: Top {}%\nDominant Frequencies: {}\n\n".format(in_dir,f_count,a_time,s_time,d_thresh,d_freq_count))
            f.write("Parameter info:\nLambda S: {}\nTau_min: {}\nAva_angle: {}\nHeight: {} Depth: {} Length: {}\nTotal Frequencies (amplitude above {}): {}\n\n".format(lam_s,tau,Ava,h,d,l,thresh,freqs))

            #Obtain trendline/function equations from analysis
            #Summary fields: "Total Time","X","Y","SD_Amplitude","Avg_Amplitude","SD_Phase_Velocity","Avg_Phase_Velocity","Wavelength","Velocity","Avg_Aspect_Ratio"
            #Top amp fields: "Time","Dominant","X","Y","Amplitude","Phase","PhaseVelocity","Wavelength","Amp/Wave"
            dispersion = str(np.poly1d(np.polyfit(summary_data["Wavelength"],summary_data["Velocity"],1)))
            fit = curve_fit(func,top_amp_data["Time"],(top_amp_data["X"]+(top_amp_data["Y"])*d))[0]
            wave_time = "({})e^(-t/({})) - ({})t + ({})".format(fit[0],fit[1],fit[2],fit[3])
            f.write("Dispersion:{}\n\nWavelength/Time:\n{}\n\n".format(dispersion,wave_time))
        
        pd.set_option('display.max_columns',None)
        pd.set_option('display.max_rows',None)
        pd.set_option('expand_frame_repr',False)
        f.write("Summary data below:\n{}\n".format(summary_data))
        f.write("\nTop Frequencies Per Timestep:\n{}\n".format(top_amp_data))
    except Exception as e:
        print(e.message)
        return
    finally:
        f.close()

    print("Summary for directory at: {}{}".format(directory,filename))

#Reads all files (simulation output) in directory with specific extension and returns it as a sorted list
#dir_path -> The directory containing data files
#pref -> The standard prefix of datafiles
#datatype -> The type of data stored in each files (int or float data)
def read_directory(dir_path,pref,par_ext,datatype,skip_files,verbose=True):

    MAX_ERRORS = 10 #Stop reading in files if max is reached
    
    if verbose:
        print("Sorting through data files...")
    #Get array list of all files to open
    files = []
    par_file_path = ""
    try:
        for f in os.listdir(dir_path):
            if f.startswith(pref):
                files.append(os.path.join(dir_path,f))
            elif f.endswith(par_ext):
                par_file_path = dir_path+"/"+f
    except:
        if verbose:
            print("An error occured, check directory path and the read/write permissions.\nDirectory: {}".format(dir_path))    
        return
    
    files.sort()
    count = 0.0
    max = len(files)/skip_files

    error_count = 0 #Track how many error files encountered
    show_errors = True
    #Create list of numpy arrays containing data for each file
    if verbose:
        print("Reading data files..")
    all_data = []
    for i, f in enumerate(files):
        if i % skip_files == 0:

            if error_count > MAX_ERRORS and verbose:
                print("Max file reading errors has been reached. {} files had errors in data.".format(error_count))
                print("Total files read successfully: {}. Skipping to analysis process...".format(i+1-error_count))
                break
            np_arr = read_data(f,datatype,show_errors)
            if np_arr is None:
                error_count += 1
                # PY3 Note: Debugging Purposes, show all errors
                # if show_errors:
                    # show_errors = False
            else:
                all_data.append(np_arr)
            
            if verbose:
                count += 1
                progress = (count / max) * 100.0
                sys.stdout.write('\r[{}] {}%'.format('#'*int(progress/5), round(progress,2)))
                sys.stdout.flush()
    if error_count > 0 and verbose:
        print("{} files caused errors when reading in data.".format(error_count))

    return [all_data, par_file_path]

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
    # PY3 Note: need indeces to be integers
    return np.fft.fft2(fft_data)[0:int(dpx/2),0:int(dpy/2)]*2/(dpx*dpy)

# Run fft2d_analyze on a list of many files
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
    top_amps = []
    max_a = np.argmax
    shape = all_fft[0].shape
    coord = np.unravel_index
    append_t = top_amps.append
    append_a = all_amps.append
    amp = np.abs
    for data in all_fft:
        a = amp(data)
        append_t(coord(max_a(a),shape))
        append_a(a)
    
    return all_amps, top_amps

#Calculates all the phases of the fft_analysis results
#all_fft -> fft results taken from all the data in a directory
def all_phases(all_fft):
    
    all_phases = []
    for data in all_fft:
        all_phases.append(np.angle(data)%(2*np.pi))
    return all_phases

#Locates the dominant frequencies within amplitude data by filtering with a threshold value
#Returns a list of dominant frequencies as well as a full list of all data of interest
#all_amps -> a list of numpy arrays containing amplitude data from fft results
#threshold -> the value used to filter out and detect dominant frequencies
def get_dominant_freqs(all_amps, threshold, noise_thresh):
    dominant_freqs = []
    for amps in all_amps:
        top = top_percent(amps, threshold)
        if len(top) == 1 and top[0][0] > noise_thresh:
            coords = (top[0][1],top[0][2])
            if coords not in dominant_freqs:
                dominant_freqs.append(coords)
    return dominant_freqs

#Creates a list of frequencies which meet the minimum amplitude threshold set, includes dominant frequencies
def purge_noise_freqs(all_amps, threshold):
    freqs = []
    l, w = all_amps[0].shape

    for amps in all_amps:
        if np.amax(amps) < threshold:
            continue
        else:
            for y in range(0,l-1):
                for x in range(0,w-1):
                    if amps[y][x] >= threshold and not (x,y) in freqs:
                        freqs.append((x,y))

    return freqs

#Get the top values based on a percentile threshold. Returns a list of values [value, x, y]
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

#Filters out all values below specified threshold, returns a list of values [value, x, y]
#values -> 2d array containing the values to filter
#threshold -> values that are below threshold are removed
def value_threshold(values, threshold):
    
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

#Calculates all of the phase velocities of all files
#all_phases -> a list of 2d numpy arrays of all phase data
#all_amps -> a list of 2d numpy arrays with all amplitudes
#time_delta -> the time change between one array and another
def get_all_velocities(all_phases, all_amps, time_delta):
    
    p_velocities = []
    
    for i, phase in enumerate(all_phases[1:]):
        diff = ((phase-all_phases[i]+2*np.pi) % 2*np.pi) - 2*np.pi
        p_velocities.append(diff/time_delta)

    return p_velocities

#Calculates features at a specific x,y frequency
#time -> the current time to calculate
#x, y -> frequencies to calculate
#amps, phases, velocities -> set of amplitudes, phases and phase velocities for that specific time
#d_freqs -> dominant frequencies that were discovered
def get_data_at_freq(time,x,y,amps,phases,velocities,d_freqs):
    dom = (x,y) in d_freqs
    flat_coord = np.ravel_multi_index((x,y),amps.shape)
    amp = amps[x][y]
    phases = phases[x][y]
    pv = velocities[x][y]
    if x > 0:
        wave = 1.0/x
    else:
        wave = 1
    aspect = amp/wave
    data = [time,flat_coord,x,y,dom,amp,phases,pv,wave,aspect]
    return data

#Creates a pandas dataframe containing all data for a specific frequency at all times of simulation
#time_step -> the amount of time per data file
#x,y -> frequency that the dataframe is about
#all_amps, all_phases etc -> lists of 2d numpy arrays containing all the data taken from files
#d_freqs -> dominant frequencies
def build_frame(time_step,x,y,all_amps,all_phases,all_velocities,d_freqs):

    stats = []
    total_time = len(all_amps)

    for t in range(1,total_time):
        stats.append(get_data_at_freq(t*time_step,x,y,all_amps[t],all_phases[t],all_velocities[t-1],d_freqs))

    return pd.DataFrame(stats,columns=['Time','Flat-Coord','X','Y','Dominant','Amplitude','Phase','PhaseVelocity','Wavelength','Amp/Wave'])

#Creates a master dataframe that contains all data frames from all frequencies which have a max amplitude equal or above the threshold
#Also creates a summary dataframe with information about the entire set
def build_all_frames(freqs,time_step,top_amps,all_amps,all_phases,all_velocities,d_freqs):
    
    frames = []
    summary_data = []
    dom_freq_stats = []
    top_amp_stats = []
    #Build frames for each x, y frequency 
    for coords in freqs:
        x, y = coords
        
        frame = build_frame(time_step,x,y,all_amps,all_phases,all_velocities,d_freqs)
        frames.append(frame)
        
        #Use frame to get some more information about this frequency
        avgPV = frame["PhaseVelocity"].mean()
        stdPV = frame["PhaseVelocity"].std()
        wave = frame.at[0,"Wavelength"]
        velocity = avgPV*wave
        avgAmp = frame["Amplitude"].mean()
        stdAmp = frame["Amplitude"].std()
        avgAspec = avgAmp / avgPV
        summary_data.append([frame.at[0,"Dominant"],frame.at[0,"X"],frame.at[0,"Y"],stdAmp,avgAmp,stdPV,avgPV,wave,velocity,avgAspec])

    total_time = len(all_amps)
    for t in range(1,total_time):
        x, y = top_amps[t]
        top_amp_stats.append(get_data_at_freq(t*time_step,x,y,all_amps[t],all_phases[t],all_velocities[t-1],d_freqs))

    
    #Create summary dataframes
    summary_frame = pd.DataFrame(summary_data,columns=["Dominant","X","Y","SD_Amplitude","Avg_Amplitude","SD_Phase_Velocity","Avg_Phase_Velocity","Wavelength","Velocity","Avg_Aspect_Ratio"])
    top_amp_frame = pd.DataFrame(top_amp_stats,columns=["Time","Flat_X-Y","X","Y","Dominant","Amplitude","Phase","PhaseVelocity","Wavelength","Amp/Wave"])

    #Concatenate all frames into single large data frame and return along with summary frame
    return [pd.concat(frames),summary_frame,top_amp_frame]

#Returns the fft results as a Pandas dataframe 
#timestep -> the current frame/time of the data being passed
#threshold -> the threshold value of data to track when providing the stats
#amps -> 2d array that contains the amplitude results of fft
#phases -> 2d array that contains the phase results of fft
#d_freqs -> a list of the dominant frequencies found among all data
def get_stats(timestep,threshold,amps,phases,d_freqs):

    #get stats for top values
    top_vals = top_values(amps, threshold, d_freqs)

    stats = []

    for i, data in enumerate(top_vals):
        amp = data[0]
        x = data[1]
        y = data[2]
        dom = data[3]
        phase = phases[y][x]
        if x > 0:
            wave = 1.0/x
        else:
            wave = 1
        aspect = amp/wave
        stats.append([timestep,dom,x,y,amp,phase,wave,aspect])
    
    return pd.DataFrame(stats,columns=['Time','Dominant Freq.','X','Y','Amplitude','Phase','Wavelength','Amp/Wave'])

#Returns one large dataframe containing all the results taken from fft analysis of all files in directory
#skip_val -> The amount of timesteps skipped.
#threshold -> the threshold value of data to track when providing the stats
#d_freqs -> a list of the dominant frequencies found among all data
#all_phases -> a list of numpy arrays containing the phase data taken from all the files
#all_amps -> a list of numpy arrays conatining all the amplitude data
def get_all_stats(skip_val, threshold, d_freqs, all_phases, all_amps):
    
    frames = []
    for i, amp in enumerate(all_amps):
        frames.append(get_stats(i*skip_val,threshold,amp,all_phases[i],d_freqs))

    #Create large dataframe, skip first one since its all noise
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

    import imageio

    if interval == 0:
        return 0

    count = 0.0
    max = len(data)/interval
    images = []
    cmin = np.amin(data)
    cmax = np.amax(data)*.9
    for i, d in enumerate(data):
        if i % interval == 0:
            count += 1
            progress = count / max * 100
            sys.stdout.write('\r[{}] {}%'.format('#'*int(progress/5), round(progress,2)))
            sys.stdout.flush()
            fig = plot_data(d,cmin,cmax,graph_type,fig_size,x_label,y_label,z_label,title.format(id=i))
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
def plot_data(data,min,max,type='wire',fig_size=(10,10),x_label='x',y_label='y',z_label='Amplitude',title='Graph of FFT2d Amplitude Spectrum'):
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
        ax = fig.add_subplot(111)
        pl.contourf(data[0:15,0:30],vmin=min,vmax=max,cmap=pl.cm.get_cmap("BuPu"))
    else:
        ax.plot_wireframe(data_x,data_y,data)

    return fig
    
#Saves the figure as a png image
def save_to_png(figure, fname):
    figure.savefig(fname+'.png',bbox_inches='tight')

#Performs fft analysis on specified directory
def analyze_directory(dir_name, output_dir, base_pref, par_ext, output_name, image_interval, skip_files, verbose=True, shared_q=None):

    MAIN_DATA_DIR = dir_name
    PARAMETER_FILE = par_ext
    SKIP_FILES = skip_files
    PNG_OUTPUT_DIR = output_dir + "png_output/"
    DATA_OUTPUT_DIR = output_dir
    CSV_OUTPUT_NAME = output_name + "_{}.csv".format(int(t.time()))
    SUMMARY_NAME = output_name + "_summary_{}.txt".format(int(t.time()))
    GIF_OUTPUT_NAME = output_name + ".gif"
    GRAPH_TYPE = 'cont' #Options available to use, 'surf'->surface, 'wire'->wireframe, 'scat'->scatter, 'cont'->contour
    XLABEL = 'x'
    YLABEL = 'y'
    ZLABEL = 'Amplitude'
    TITLE = 'Graph of FFT2d Amplitudes at t={id}' #id is replaced by index of snapshot
    FIG_SIZE = (10,10)
    SNAPSHOT_INTERVAL = image_interval #A png image of graph is saved at each interval
    BASE_FILE_NAME = "ALTI{:05d}_t0"
    BASE_PREF = base_pref
    THRESHOLD = 0.8
    AMP_THRESHOLD = 0.7 #Should be around 1
    DATA_TYPE = int
    TIME_DELTA = 10 #Time change between data files
    
    try:
        #Check parent directories exists
        directories = [MAIN_DATA_DIR]
        for d in directories:
            if not os.path.exists(d):
                err = "The specified path: {} was not found. Analysis cancelled.".format(d)
                print(err)
                if shared_q is not None:
                    shared_q.put({"Error":err,"Summary":SUMMARY_NAME})
                    shared_q.close()
                return 1

        #Make directories if needed
        if image_interval > 0 and not os.path.exists(PNG_OUTPUT_DIR):
            # PY3 Note: 0777 permissions are default
            os.makedirs(PNG_OUTPUT_DIR)
        
        if not os.path.exists(DATA_OUTPUT_DIR):
            # PY3 Note: 0777 permissions are default            
            os.makedirs(DATA_OUTPUT_DIR)

        #Get all data from the main directory
        t0 = t.time()
        all_data, par_file_path = read_directory(MAIN_DATA_DIR,BASE_PREF,par_ext,DATA_TYPE,SKIP_FILES,verbose)
        t1 = t.time()

        #Check input files exists
        file_count = len(all_data)
        if file_count == 0:
            err = "No files with the prefix: '{}' were found. Analysis cancelled.".format(BASE_PREF)
            if verbose:
                print(err)
            if shared_q is not None:
                shared_q.put({"Error":err,"Summary":SUMMARY_NAME})
                shared_q.close()
            return 1
        t_read = t1-t0
        
        if verbose:
            print("\rTotal files read: {} Read time: {}s\nCalculating fft2d data on all files...".format(file_count,t_read))
        
        #Get all fft2d data for calculating
        t0 = t.time()
        all_fft2d = all_fft2d_analysis(all_data)
        t1 = t.time()
        t_fft2d = t1-t0
        if verbose:
            print("FFT calculations time: {}s\nCalculating phase data...".format(t_fft2d))

        #Get all phase data
        t0 = t.time()
        all_phase_data = all_phases(all_fft2d)
        t1 = t.time()
        t_phases = t1 - t0
        if verbose:
            print("Phases calculation time: {}s\nCalculating amplitudes...".format(t_phases))

        #Get all amplitude data
        t0 = t.time()
        all_amps, top_amps = all_amplitudes(all_fft2d)
        t1 = t.time()
        t_amps = t1-t0
        if verbose:
            print("Amplitude calculation time: {}s\nCalculating phase velocities...".format(t_amps))
            
        #Get all phase velocities
        t0 = t.time()
        all_velocities = get_all_velocities(all_phase_data, all_amps, TIME_DELTA)
        t1 = t.time()
        t_velocities = t1-t0
        if verbose:
            print("Phase Velocity calculation time: {}s\nFinding dominant frequencies...".format(t_velocities))

        #Find dominant frequencies
        t0 = t.time()
        d_freqs = get_dominant_freqs(all_amps,THRESHOLD,AMP_THRESHOLD)
        t1 = t.time()
        t_freqs = t1-t0
        if verbose:
            print("{} dominant frequencies found at threshold {}%, time: {}s\nDeriving all data and writing to csv...".format(len(d_freqs),THRESHOLD*100,t_freqs))
        
        #Concatenate and save data results
        t0 = t.time()
        freqs = purge_noise_freqs(all_amps,AMP_THRESHOLD)
        master_frame, summary, top_amp_summary = build_all_frames(freqs,TIME_DELTA,top_amps,all_amps,all_phase_data,all_velocities,d_freqs)
        master_frame.to_csv(DATA_OUTPUT_DIR + CSV_OUTPUT_NAME)
        analysis_time = t_read + t_fft2d + t_phases + t_amps + t_freqs

        #Change permissions to read/write for all and directories
        err = "Error changing CSV output permissions. File: {}".format(DATA_OUTPUT_DIR + CSV_OUTPUT_NAME)
        try:
            os.chmod(DATA_OUTPUT_DIR + CSV_OUTPUT_NAME, 0o777)
            err = "Issue writing summary file"
            run_stats = [MAIN_DATA_DIR, file_count, analysis_time, len(d_freqs),THRESHOLD,AMP_THRESHOLD]
            write_summary(run_stats,DATA_OUTPUT_DIR,SUMMARY_NAME,par_file_path,summary,top_amp_summary)
            err = "Error changing summary file permissions."
            os.chmod(DATA_OUTPUT_DIR + SUMMARY_NAME, 0o777)
        except:
            if verbose:
                print(err)
            if shared_q is not None:
                shared_q.put({"Error":err,"Summary":SUMMARY_NAME})
                shared_q.close()
            return 1

        t1 = t.time()
        t_stats = t1 - t0
                    
        if verbose:
            print("Data calculation and concatenation time: {}".format(t_stats))

        #Graph resulting data at specific intervals and save as images to directory, create GIF of pngs
        if image_interval > 0 and can_plot:

            if verbose:
                print("Plotting data at intervals of {} and creating PNG images...".format(t_stats,SNAPSHOT_INTERVAL))
            t0 = t.time()
            im_count = graph_all(SNAPSHOT_INTERVAL,PNG_OUTPUT_DIR,DATA_OUTPUT_DIR+GIF_OUTPUT_NAME,BASE_FILE_NAME,all_amps,GRAPH_TYPE,FIG_SIZE,XLABEL,YLABEL,ZLABEL,TITLE)
            t1 = t.time()
            t_plot = t1-t0
            if verbose:
                print("\r{} PNG's created in: {}s.\nGIF animation complete.\nAnalysis process complete!\nTotal time: {}s".format(im_count,t_plot,analysis_time + t_stats + t_plot))
        elif verbose:
            print("Analysis of direcory: {} complete! Analysis time: {}s".format(MAIN_DATA_DIR,analysis_time + t_stats))

        if shared_q is not None:
            shared_q.put({"Error":None,"Summary":summary})
            shared_q.close()
    except Exception as e:
        err = e.message
        if verbose:
            print(err)
        if shared_q is not None:
            shared_q.put({"Error":err,"Summary":SUMMARY_NAME})
            shared_q.close()

def analyze_many_dir(main_dir, output_dir, base_pref, par_ext, img_int, skip_files):

    t0 = t.time()
    
    #Get list of sub directories to analyze
    dirs = []
    try:
        for d in os.listdir(main_dir):
            if os.path.isdir(main_dir+d):
                dirs.append(d)
    except:
        print("An error occured, check that main directory path is correct.")

    proc = []
    shared_q = mp.Queue()
    for d in dirs:
        main_d = main_dir + d
        print("Processing directory: {}".format(main_d))
        p = mp.Process(target=analyze_directory,args=(main_d,output_dir,base_pref,par_ext,d,img_int,skip_files,False,shared_q))
        p.start()
        proc.append(p)

    for p in proc:
        p.join()

    t1 = t.time()
    if len(dirs) >= 1:
        print("All directories analyzed! Total Time: {}".format(t1-t0))

        #Create a summary of summaries for all successful processes
        try:
            f = open(output_dir+"MAIN_SUMMARY.TXT","a+")
            while not shared_q.empty():
                p_data = shared_q.get()
                if p_data["Error"] is None:
                    f.write(p_data["Summary"]+"\n")
                else:
                    print("Creating {} file not successful.\nError: {}".format(p_data["Summary"],p_data["Error"]))
            f.close()
            os.chmod(output_dir+"MAIN_SUMMARY.TXT", 0o777)
        except:
            print("An error occured while writing main summary file.")
    else:
        print("No directories were analyzed. Check the main directory path.\nMain directory used: {}".format(main_dir))

# KK: how is this different from analyse directory?
# Should it exist as an independent function?
def plot_only(dir_name, output_dir, base_pref, par_ext, output_name, image_interval, skip_files, verbose=True):
    MAIN_DATA_DIR = dir_name
    PARAMETER_FILE = par_ext
    SKIP_FILES = skip_files
    PNG_OUTPUT_DIR = output_dir + "png_output/"
    DATA_OUTPUT_DIR = output_dir
    CSV_OUTPUT_NAME = output_name + ".csv"
    SUMMARY_NAME = output_name + "_summary.txt"
    GIF_OUTPUT_NAME = output_name + ".gif"
    GRAPH_TYPE = 'cont' #Options available to use, 'surf'->surface, 'wire'->wireframe, 'scat'->scatter, 'cont'->contour
    XLABEL = 'x'
    YLABEL = 'y'
    ZLABEL = 'Amplitude'
    TITLE = 'Graph of FFT2d Amplitudes at t={id}' #id is replaced by index of snapshot
    FIG_SIZE = (10,10)
    SNAPSHOT_INTERVAL = image_interval #A png image of graph is saved at each interval
    BASE_FILE_NAME = "ALTI{:05d}_t0"
    BASE_PREF = base_pref
    THRESHOLD = 0
    AMP_THRESHOLD = 0 #Should be
    DATA_TYPE = int
    TIME_DELTA = 10 #Time change between data files
    
    #Check parent directories exists
    directories = [MAIN_DATA_DIR]
    for d in directories:
        if not os.path.exists(d):
            print("The specified path: {} was not found. Analysis cancelled.".format(d))
            return 1

    #Make directories if needed
    if image_interval > 0 and not os.path.exists(PNG_OUTPUT_DIR):
        # PY3 Note: 0777 permissions are default
        os.makedirs(PNG_OUTPUT_DIR)
    
    if not os.path.exists(DATA_OUTPUT_DIR):
        # PY3 Note: 0777 permissions are default
        os.makedirs(DATA_OUTPUT_DIR)

    #Get all data from the main directory
    t0 = t.time()
    all_data, par_file_path = read_directory(MAIN_DATA_DIR,BASE_PREF,par_ext,DATA_TYPE,SKIP_FILES,verbose)
    t1 = t.time()

    #Check input files exists
    file_count = len(all_data)
    if file_count == 0:
        print("No files with the prefix: '{}' were found. Analysis cancelled.".format(BASE_PREF))
        return 1
        t_read = t1-t0
        
        if verbose:
            print("\rTotal files read: {} Read time: {}s\nCalculating fft2d data on all files...".format(file_count,t_read))

    #Get all fft2d data for calculating
    t0 = t.time()
    all_fft2d = all_fft2d_analysis(all_data)
    t1 = t.time()
    t_fft2d = t1-t0
    if verbose:
        print("FFT calculations time: {}s\nCalculating amplitude data...".format(t_fft2d))

    #Get all amplitude data
    t0 = t.time()
    all_amps = all_amplitudes(all_fft2d)
    t1 = t.time()
    t_amps = t1-t0
    if verbose:
        print("Amplitude calculation time: {}s".format(t_amps))

    #Graph resulting data at specific intervals and save as images to directory, create GIF of pngs
    if image_interval > 0 and can_plot:
        if verbose:
            print("Graphing data")
        im_count = graph_all(SNAPSHOT_INTERVAL,PNG_OUTPUT_DIR,DATA_OUTPUT_DIR+GIF_OUTPUT_NAME,BASE_FILE_NAME,all_amps,GRAPH_TYPE,FIG_SIZE,XLABEL,YLABEL,ZLABEL,TITLE)

    print("\nPlotting done! {} images graphed.".format(im_count))

#directory -> The directory that holds the log files to analyze
#image_interval -> A graph will be made and saved at every interval. E.g 50 = every 50th data file will be graphed.
#Note: if image interval is set to 0, no graphs are made.
def main(directory="input_data/ALT_DATA1/",output_dir="ALT_DATA1_OUT",filename="",image_interval=100, verbose=False, plot=False):

    if verbose in ["True", "TRUE", "true", "1", "t"]:
        verbose = True
    else:
        verbose = False

    if plot in ["True", "TRUE", "true", "1", "t"]:
        plot = True
    else:
        plot = False
        
    if plot:
        plot_only(directory,output_dir,"ALTI",".par",filename,image_interval,1,True)
    else:
        if filename=="":
            analyze_many_dir(directory,output_dir,"ALTI",".par",0,1)
        else:
            analyze_directory(directory,output_dir,"ALTI",".par",filename,image_interval,1,verbose)
    
args = sys.argv
argcount = len(args)

if argcount > 6:
    main(args[1],args[2],args[3],int(args[4]),args[5],args[6])
elif argcount > 5:
    main(args[1],args[2],args[3],int(args[4]),args[5])
else:
    print("Not enough arguments passed to function.")

# PY# Note: Testing Calls
# main(directory="input_data/ALT_DATA1/",output_dir="ALT_DATA1_OUT",filename="ALT_DATA1",image_interval=1, verbose=True, plot=True)

# main(directory="input_data/Series", output_dir="ALT_DATA1_OUT",filename="Series",image_interval=1, verbose=True, plot=True)
