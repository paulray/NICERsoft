"""
A script to fill data gaps in a way that retains the noise characteristic of the live data set and preserves AVAR curves. 
An implementation of the algorithm described in:
Howe and Schlossberger, "Characterizing Frequency Stability Measurements having Multiple Data Gaps"
"""

import sys
import numpy as np
import copy
import csv
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from itertools import groupby
from operator import itemgetter

def fillgaps(datafile, savefile = "NaN", plotresults = True):
    """A method to fill gaps in data and preserve data noise characteristics.

    Parameters:
    datafile: (string) a file path ending in ".csv" to a comma delimited csv file with two columns formatted as 
              <time,equispaced>, <data>

    savefile: (string) a file path ending in ".csv" to which to save the imputed results
    
    plotresults: (Bool) whether or not to make a plot of the input data and filled data

    """

    y = []
    x = np.array([])
    if datafile[-4:] == ".csv":  # if .csv file entered
        with open(datafile, "r") as f:
            lines = f.readline()
            while lines:
                vals = lines.split(",")
                x = np.append(x, float(vals[0]))
                try:
                    y.append(float(vals[1]))
                except:
                    y.append(np.nan)
                lines = f.readline()
    elif datafile[-4:] == ".txt":  # if .txt file entered
        with open(datafile, "r") as f:
            lines = f.readline()
            while lines:
                vals = lines.strip().split()
                if vals != []:  # if not an empty line
                    x = np.append(x, float(vals[0]))
                    y.append(float(vals[1]))
                lines = f.readline()

    # create equispaced array for x-axis
    if len(x) >= 2:
        diff = x[1] - x[0]
        for i in range(len(x) - 1):
            if x[i + 1] - x[i] < diff and x[i + 1] - x[i] > 0:
                diff = x[i + 1] - x[i]
        numintervals = int((x[-1] - x[0]) / diff)
        xindices = np.array([0], dtype=int)
        for i in range(1, len(x)):
            for j in range(1, numintervals + 1):
                if x[i] - (x[0] + diff * j) < diff:
                    xindices = np.append(xindices, j)
                    break
    data = np.zeros(xindices.max() + 1) * np.nan
    data[xindices] = y

    xfilled = np.linspace(x[0], x[-1], num=int(numintervals + 1))
    if plotresults:
        fig, ax = plt.subplots(2, figsize=(12, 8))
        ax[0].set_title("Original Data")
        ax[0].scatter(x, y, marker = '.')
        ax[0].grid(True)
   

    # find indeces of data gaps and largest continuous run of data
    gap_indeces = np.where(np.isnan(data))[0]
    # check if the above numpy.ndarray empty - print message and exit
    if gap_indeces.size == 0:
        print("The data has no gaps. Exiting...")
        return

    gap_total = 1
    start_gap_inx = gap_indeces[0]
    maxdiff = start_gap_inx  # run before first gap
    gap_to_impute = 1
    # initial value
    gaps = {
        1: (start_gap_inx, start_gap_inx)
    }  # {gap #: (starting index, ending index)}
    for i in range(len(gap_indeces) - 1):
        diff = gap_indeces[i + 1] - gap_indeces[i]
        if diff > maxdiff:
            maxdiff = diff
            gap_to_impute = (
                gap_total + 1
            )  # gap to the right of the largest run will be the first to impute
        if diff > 1:  # new gap reached
            gaps[gap_total] = (start_gap_inx, gap_indeces[i])
            start_gap_inx = gap_indeces[i + 1]
            gap_total = gap_total + 1
            # check if new gap also at last index (1-index gap)
            if i == len(gap_indeces) - 2:
                gaps[gap_total] = (start_gap_inx, start_gap_inx)
        if i == len(gap_indeces) - 2:  # only one jump in data set
            gaps[gap_total] = (start_gap_inx, gap_indeces[i + 1])

    # impute single-point gaps
    i = 1
    while i <= gap_total:
        # for i in range(1, gap_total):
        if gaps[i][0] == gaps[i][1]:
            # if data on both sides of gap, take average of previous and following data point
            if gaps[i][0] != 0 and gaps[i][0] != len(data) - 1:
                data[gaps[i][0]] = (data[gaps[i][0] - 1] + data[gaps[i][0] + 1]) / 2
            elif gaps[i][0] == 0:
                data[0] = data[1]
            else:
                data[-1] = data[-2]
            for k in range(i, gap_total):
                gaps[k] = gaps.pop(k + 1)
            gap_total = gap_total - 1
        else:
            i = i + 1

    # check if first gap_to_impute now beyond the range of total gaps
    if gap_to_impute > gap_total:
        gap_to_impute = 1  # set to first gap - if only 1 gap, automatic choice, otherwise starting point
        if gap_total > 1:
            maxdiff = gaps[1][0]  # run before first gap
            for i in range(1, gap_total):
                diff = gaps[i + 1][0] - gaps[i][1]
                if diff > maxdiff:
                    maxdiff = diff
                    gap_to_impute = i + 1

    # impute gaps to the right of this run, then reverse and do left 
    # i. Type 2 - reflect+invert

    while gap_total != 0:
        while gap_to_impute <= gap_total and gap_to_impute != 0:
            gap_to_impute, gap_total = fill(
            data, gaps, gap_total, gap_to_impute, xfilled 
            )
            if gap_to_impute is None:
                return

        # spacially reverse dataset, continue until beginning reached
        for i in range(gap_to_impute, 0, -1):
            gap_to_impute, gap_total = fill(
                data, gaps, gap_total, i, xfilled, reverse=True
            )
            if gap_to_impute is None:
                return

        # reset gap_to_impute if still gaps left
        maxdiff = 0
        for i in range(1, len(gaps)):
            diff = gaps[i + 1][0] - gaps[i][1]
            if diff > maxdiff:
                gap_to_impute = (
                    i + 1
                )  # gap to the right of the largest continuous run
    if savefile != "NaN":
        with open(savefile, "w",  newline='') as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(zip(xfilled, data))  # write results into csv file in same format as input
    if plotresults:
        ax[1].scatter(xfilled, data, marker='.')
        ax[1].grid(True)
        ax[1].set_title("Recovered Data")
        plt.show()  # displays graphs of data before/after imputation



def filterfunction(data):
    """
    Takes a continuous block of data and multiplies it in the frequency domain by exp(-2^3 f/T)
    (The fft shows the positive frequencies followed by flipped negative frequencies, so we flip the exponential filter and append itself to its end)
    """
    #pad dataset
    data = np.append(np.flip(data),np.append(data,np.flip(data)))
    #fft
    fourierdata = fft(data.T).T
    #multiply by filter
    length = np.size(data)
    if (length%2)==0:
        filterform = np.exp(-np.arange(0,length/2)*(16/length))
        filterform = np.append(filterform,np.flip(filterform))
    else:
        filterform = np.exp(-np.arange(0,(length-1)/2)*(16/length))
        filterform = np.append(filterform,np.append(filterform[-1],np.flip(filterform)))
    filteredfourier = np.multiply(fourierdata.reshape(length,1),filterform.reshape(length,1))
    #ifft
    filtereddata = np.real(ifft(filteredfourier.T).T)
    #unpad
    filtereddata = filtereddata[int(length/3):int(2*length/3)]
    return filtereddata
def filter_all(data):
    """
    Takes a set of data with gaps and applies the filter to each continuous block of data within the dataset
    """
    #filter individual blocks of continuous data
    datablocks = [list(v) for k,v in groupby(data,np.isfinite) if k]
    filtereddatablocks = [filterfunction(v) for v in datablocks]
    if np.any(np.isnan(data)): #if there are gaps
        #get information about gaps
        gap_indeces = np.where(np.isnan(data))[0]
        gap_locs = []
        gap_lengths = []
        for k, g in groupby(enumerate(gap_indeces), lambda ix:ix[1]-ix[0]):
            temprun = list(map(itemgetter(1), g))
            gap_lengths.append(np.size(temprun))
            gap_locs.append(temprun[0])
        #build filtered data with gaps
        filteredrun = np.array([]);
        if gap_locs[0]==0: #if the first space is a gap, append and remove
            filteredrun = np.append(filteredrun, np.full(gap_lengths[0],np.NaN) )
            gap_locs = np.delete(gap_locs,0)
            gap_lengths = np.delete(gap_lengths,0)
        for i in range(np.size(filtereddatablocks)):
            filteredrun = np.append(filteredrun,filtereddatablocks[i])
            if i<np.size(gap_locs):
                filteredrun = np.append(filteredrun,np.full(gap_lengths[i],np.NaN))
    else: #if no gaps
       filteredrun = filterfunction(data)
    return filteredrun

def check_right(data, gaps, curgap, curgap_size, curgap_num, gap_total):
    """A method to check the right-hand side of the data to obtain the correct number of points to fill the gap.

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            curgap - tuple, (start, end) indeces for current gap
            curgap_size - int, size of current gap
            curgap_num - int, number of current gap
            gap_total - int, total number of gaps
    """
    if (
        len(data) - 1 - curgap[1]
    ) >= curgap_size:  # if enough data to the right of the gap
        if curgap_num == gap_total:  # if on the last gap
            pts_to_reflect = copy.deepcopy(
                data[(curgap[1] + 1) : (curgap[1] + 1 + curgap_size)]
            )  # take the next number of data points required to fill up the gap
            return pts_to_reflect, "right"

        elif (gaps[curgap_num + 1][0] - 1) - (
            curgap[1] + 1
        ) >= curgap_size:  # if continuous run between this gap and the next has enough points to fill the gap
            pts_to_reflect = copy.deepcopy(
                data[(curgap[1] + 1) : (curgap[1] + 1 + curgap_size)]
            )  # take the next number of data points required to fill up the gap
            return pts_to_reflect, "right"

        else:
            # not enough data points between this gap and the next, either print error message and exit or impute from left and right sides
            pts_to_reflect = check_both_sides(
                data, gaps, curgap, curgap_size, curgap_num, gap_total
            )
            if pts_to_reflect is None:
                return None, None
            else:
                return pts_to_reflect, "both"
    else:
        pts_to_reflect = check_both_sides(
            data, gaps, curgap, curgap_size, curgap_num, gap_total
        )
        if pts_to_reflect is None:
            return None, None
        else:
            return pts_to_reflect, "both"


def check_both_sides(data, gaps, curgap, curgap_size, curgap_num, gap_total):
    """A method to check both sides of the data to obtain the correct number of points to fill the gap (half left, half right).

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            curgap - tuple, (start, end) indeces for current gap
            curgap_size - int, size of current gap
            curgap_num - int, number of current gap
            gap_total - int, total number of gaps
    """
    # either print error message and exit or impute from left and right sides
    if curgap_num == 1:  # if this is the first gap in the dataset
        if (
            (gap_total == 1)
            and (curgap[0] >= curgap_size / 2)
            and (curgap[1] <= len(data) - (curgap_size / 2))
        ):  # if only gap
            pts_to_reflect_left = copy.deepcopy(
                data[int(curgap[0] - curgap_size / 2) : (curgap[0])]
            )
            pts_to_reflect_right = copy.deepcopy(
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) + 1]
            )
            return pts_to_reflect_left, pts_to_reflect_right
        elif (
            (gap_total != 1)
            and (curgap[0] >= curgap_size / 2)
            and ((gaps[curgap_num + 1][0] - 1) - (curgap[1] + 1) >= curgap_size / 2)
        ):  # if enough space between this gap and the next
            pts_to_reflect_left = copy.deepcopy(
                data[int(curgap[0] - curgap_size / 2) : (curgap[0])]
            )
            pts_to_reflect_right = copy.deepcopy(
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) ]#+ 1]
            )
            return pts_to_reflect_left, pts_to_reflect_right
        else:  # not enough space on one/both sides
            if gap_total == 1:
                print("Unable to fill all gaps, not enough data. Exiting...")
                return "Done"
            else:
                return None
    elif (
        curgap[0] >= gaps[curgap_num - 1][1] + curgap_size / 2
    ):  # if enough space between previous gap and current gap
        if (curgap[1] <= len(data) - (curgap_size / 2)) and (
            curgap_num == gap_total
        ):  # if last gap and enough space to the right
            pts_to_reflect_left = copy.deepcopy(
                data[int(curgap[0] - curgap_size / 2) : (curgap[0])]
            )
            pts_to_reflect_right = copy.deepcopy(
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) ]#+ 1]
            )
            return pts_to_reflect_left, pts_to_reflect_right
        elif (curgap_num != gap_total) and (
            (gaps[curgap_num + 1][0] - 1) - (curgap[1] + 1) >= curgap_size / 2
        ):  # if enough space before next gap
            pts_to_reflect_left = copy.deepcopy(
                data[int(curgap[0] - curgap_size / 2) : (curgap[0])]
            )
            pts_to_reflect_right = copy.deepcopy(
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) ]#+ 1]
            )
            return pts_to_reflect_left, pts_to_reflect_right
        else:  # not enough space on one/both sides
            if gap_total == 1:
                print("Unable to fill all gaps, not enough data. Exiting...")
                return "Done"
            else:
                return None
    else:  # not enough space to the left of the gap
        if gap_total == 1:
            print("Unable to fill all gaps, not enough data. Exiting...")
            return "Done"
        else:
            return None


def check_left(data, gaps, curgap, curgap_size, curgap_num, gap_total):
    """A method to check the left-hand side of the data to obtain the correct number of points to fill the gap.

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            curgap - tuple, (start, end) indeces for current gap
            curgap_size - int, size of current gap
            curgap_num - int, number of current gap
            gap_total - int, total number of gaps
    """
    if curgap_num == 1:  # if this is the first gap
        pts_to_reflect = copy.deepcopy(data[(curgap[0] - curgap_size) : (curgap[0])])
        return pts_to_reflect, "left"

    elif (curgap[0] - 1) - (
        gaps[curgap_num - 1][1] + 1
    ) >= curgap_size:  # if data run between previous gap and current gap large enough
        pts_to_reflect = copy.deepcopy(data[(curgap[0] - curgap_size) : (curgap[0])])
        return pts_to_reflect, "left"

    else:  # if data run between previous gap and current gap not large enough, check data set to the right
        pts_to_reflect = check_right(
            data, gaps, curgap, curgap_size, curgap_num, gap_total
        )
        return pts_to_reflect


def fill(data, gaps, gap_total, gap_num, xcoords, reverse = False):
    """A method to fill the gap_num-th gap in the data by reflecting and inverting data on one/both sides of the gap.

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            gap_total - int, total number of gaps
            gap_num - int, number of current gap
            xcoords - np.array, the x-coordinates corresponding to the data and gaps
    """
    # find size of gap
    curgap = gaps[gap_num]
    size = (curgap[1] + 1) - curgap[0]
    side = None
    # check if there is enough data previous to this gap
    
    if reverse:
        if curgap[1] > len(data) - size:
            # if not enough previous data (right-hand since flipped), check data following gap
            pts_to_reflect, side = check_left(
                data, gaps, curgap, size, gap_num, gap_total
            )
            if pts_to_reflect is None:
                # decrement gap_num - will come back to this gap once rest of left-sided gaps filled
                gap_num = gap_num - 1
            elif type(pts_to_reflect) == tuple:
                left_pts_to_reflect = pts_to_reflect[0]
                right_pts_to_reflect = pts_to_reflect[1]
            elif pts_to_reflect == "Done":
                return None

        else:
            pts_to_reflect, side = check_right(
                data, gaps, curgap, size, gap_num, gap_total
            )
            if pts_to_reflect is None:
                # decrement gap_num - will come back to this gap once rest of left-sided gaps filled
                gap_num = gap_num - 1
            elif type(pts_to_reflect) == tuple:
                left_pts_to_reflect = pts_to_reflect[0]
                right_pts_to_reflect = pts_to_reflect[1]
            elif pts_to_reflect == "Done":
                return None

    else:
        if curgap[0] < size:
            # if not enough previous data, check data following gap
            pts_to_reflect, side = check_right(
                data, gaps, curgap, size, gap_num, gap_total
            )
            if pts_to_reflect is None:
                # increment gap_to_impute - will come back to this gap once rest of right-sided gaps filled
                gap_num = gap_num + 1
            elif type(pts_to_reflect) == tuple:
                left_pts_to_reflect = pts_to_reflect[0]
                right_pts_to_reflect = pts_to_reflect[1]

        else:
            pts_to_reflect, side = check_left(
                data, gaps, curgap, size, gap_num, gap_total
            )
            if pts_to_reflect is None:
                # increment gap_to_impute - will come back to this gap once rest of right-sided gaps filled
                gap_num = gap_num + 1
            elif type(pts_to_reflect) == tuple:
                left_pts_to_reflect = pts_to_reflect[0]
                right_pts_to_reflect = pts_to_reflect[1]
        
    # impute the gap
    if pts_to_reflect is not None:
        if type(pts_to_reflect) == tuple:
            left_pts_to_reflect = - np.flip(left_pts_to_reflect)
            filt_left_pts_to_refect = filterfunction(left_pts_to_reflect)
            right_pts_to_reflect = - np.flip(right_pts_to_reflect)
            filt_right_pts_to_reflect = filterfunction(left_pts_to_reflect)
            right_pts_to_reflect = right_pts_to_reflect + (filt_left_pts_to_refect[-1]-filt_right_pts_to_reflect[0]) #match the two reflected data parts in the middle
            if (size / 2) % 2 == 0.5:
                mid_val = (left_pts_to_reflect[0] + right_pts_to_reflect[-1]) / 2
                left_pts_to_reflect = np.append(left_pts_to_reflect,mid_val)
            pts_to_reflect = np.append(left_pts_to_reflect,right_pts_to_reflect) #combine left and right reflections
        
        else:
            # invert fill points
            pts_to_reflect = -pts_to_reflect
            # reflect the points across the gap
            pts_to_reflect = np.flip(pts_to_reflect)
        
        # Match endpoints
        filt_pts_to_reflect = filterfunction(pts_to_reflect)
        filt_data = filter_all(data)
        start_diff = filt_data[curgap[0]-1]-filt_pts_to_reflect[0]
        end_diff = filt_data[curgap[-1]+1]-filt_pts_to_reflect[-1]
        line_to_add = np.arange(start_diff,end_diff,(end_diff-start_diff)/(size-.1))
        
            
        # Insert imputed points
        data[curgap[0] : (curgap[1] + 1)] = pts_to_reflect + line_to_add
        
        # remove gap
        if gap_num == gap_total:
            gaps.pop(gap_num)
            gap_num = gap_num - 1
        else:
            for k in range(gap_num, gap_total):
                gaps[k] = gaps.pop(k + 1)
        gap_total = gap_total - 1
    return gap_num, gap_total


if __name__ == "__main__":
    fillgaps(sys.argv[1], sys.argv[2])
