"""
A script to fill data gaps in a way that retains the noise characteristic of the live data set and preserves AVAR curves. An implementation of the algorithm described in:
Howe and Schlossberger, "Strategy for Characterizing Frequency Stability Measurements having Multiple Dead Times"
"""

import sys
import argparse
import numpy as np
import copy
from scipy.stats import median_abs_deviation
import csv
import matplotlib.pyplot as plt


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
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) + 1]
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
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) + 1]
            )
            return pts_to_reflect_left, pts_to_reflect_right
        elif (curgap_num != gap_total) and (
            (gaps[curgap_num + 1][0] - 1) - (curgap[1] + 1) >= curgap_size / 2
        ):  # if enough space before next gap
            pts_to_reflect_left = copy.deepcopy(
                data[int(curgap[0] - curgap_size / 2) : (curgap[0])]
            )
            pts_to_reflect_right = copy.deepcopy(
                data[(curgap[1] + 1) : int(curgap[1] + 1 + curgap_size / 2) + 1]
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


def fill(data, gaps, gap_total, gap_num, xcoords, reverse=False):
    """A method to fill the gap_num-th gap in the data by reflecting and inverting data on one/both sides of the gap.

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            gap_total - int, total number of gaps
            gap_num - int, number of current gap
            xcoords - np.array, the x-coordinates corresponding to the data and gaps
            reverse - bool, indicated whether dataset forward-oriented (gap being filled to the right of longest continuous data run) or reversed (left of longest run)
    """
    # find size of gap
    curgap = gaps[gap_num]
    size = (curgap[1] + 1) - curgap[0]
    gap_xcoords = xcoords[curgap[0] - 1 : (curgap[1] + 2)]
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
            if len(left_pts_to_reflect) != 1:
                for j in range(0, len(left_pts_to_reflect)):
                    left_pts_to_reflect[j] = (
                        (left_pts_to_reflect[-1] + left_pts_to_reflect[-2]) / 2
                        + left_pts_to_reflect[-1]
                        - left_pts_to_reflect[j]
                    )  # add the difference between last value and current value to last value to invert
            if len(right_pts_to_reflect) != 1:
                for j in range(0, len(right_pts_to_reflect)):
                    right_pts_to_reflect[j] = (
                        right_pts_to_reflect[0] + right_pts_to_reflect[1]
                    ) / 2 + (right_pts_to_reflect[0] - right_pts_to_reflect[j])

            # reflect data
            data[curgap[0] : (curgap[0] + 1 + int(size / 2))] = left_pts_to_reflect[
                ::-1
            ]
            data[(curgap[1] - int(size / 2)) : curgap[1] + 1] = right_pts_to_reflect[
                ::-1
            ]

            # calculate slope of imputed data
            mfillpts = -(pts_to_reflect[-1] - pts_to_reflect[0]) / (
                xcoords[curgap[0] + int(size / 2)] - xcoords[curgap[0] - int(size / 2)]
            )

            # calculate slope between sides of gap
            mdesired, b = calculate_slope(
                data, curgap, gap_num, gaps, gap_total, len(pts_to_reflect), xcoords
            )

            # calculate slope to apply to imputed data
            mdiff = mdesired - mfillpts

            # create array of even x-coordinate intervals for gap
            xdiff = np.array(gap_xcoords - gap_xcoords[0])

            # apply linear slope to points
            for j in range(0, size):
                data[curgap[0] + j] = mdiff * xdiff[j] + data[curgap[0] + j] + b

            # if odd number of spaces in gap, average first left and last right values and put in middle gap value
            if (size / 2) % 2 == 0.5:
                mid_val = (left_pts_to_reflect[0] + right_pts_to_reflect[-1]) / 2
                data[curgap[1] - int(size / 2)] = mid_val

            # remove gaps
            if gap_num == gap_total:
                gaps.pop(gap_num)
                gap_num = gap_num - 1
            for k in range(gap_num, gap_total):
                gaps[k] = gaps.pop(k + 1)
            gap_total = gap_total - 1
        else:
            # invert fill points
            if side == "right":
                if len(pts_to_reflect) != 1:
                    for j in range(0, len(pts_to_reflect)):
                        pts_to_reflect[j] = (
                            pts_to_reflect[0] + pts_to_reflect[1]
                        ) / 2 + (pts_to_reflect[0] - pts_to_reflect[j])
            else:
                if len(pts_to_reflect) != 1:
                    for j in range(0, len(pts_to_reflect)):
                        pts_to_reflect[j] = (
                            (pts_to_reflect[-1] + pts_to_reflect[-2]) / 2
                            + pts_to_reflect[-1]
                            - pts_to_reflect[j]
                        )

            # reflect the points across the gap
            data[curgap[0] : (curgap[1] + 1)] = pts_to_reflect[::-1]
            # determine the slope of the reflected points
            if side == "left":
                mfillpts = -(pts_to_reflect[-1] - pts_to_reflect[0]) / (
                    xcoords[curgap[0]] - xcoords[curgap[0] - size]
                )
            else:
                mfillpts = -(pts_to_reflect[-1] - pts_to_reflect[0]) / (
                    xcoords[curgap[1] + 1 + size] - xcoords[curgap[1] + 1]
                )
            # calculate the slope between either sides of the gap
            mdesired, b = calculate_slope(
                data, curgap, gap_num, gaps, gap_total, len(pts_to_reflect), xcoords
            )
            # calculate the slope difference to determine the slope to apply to the reflected points
            mdiff = mdesired - mfillpts

            # create array of even x-coordinate intervals for the gap
            xdiff = np.array(gap_xcoords - gap_xcoords[0])

            # add linear slope to points
            for j in range(0, size):
                data[curgap[0] + j] = mdiff * xdiff[j] + data[curgap[0] + j] + b

            # remove gap
            if gap_num == gap_total:
                gaps.pop(gap_num)
                gap_num = gap_num - 1
            else:
                for k in range(gap_num, gap_total):
                    gaps[k] = gaps.pop(k + 1)
            gap_total = gap_total - 1
    return gap_num, gap_total


def calculate_slope(data, curgap, curgap_num, gaps, gap_total, fill_pts_len, xcoords):
    """A method to calculate the linear slope between both sides of a gap.

    Parameters:
            data - np.ndarray, the data being imputed
            curgap - tuple, (start, end) indeces for current gap
            curgap_num - int, number of current gap
            gaps - dict, stores the gap indeces with their respective number
            gap_total - int, total number of gaps
            fill_pts_len - int, the length of the data to fill the gap
            xcoords - np.array, the x-coordinates corresponding to the data and gaps
    """
    m = 0
    b = 0
    if fill_pts_len > 1:
        # check for outliers at the boundary and skip over as necessary
        start, end = check_boundaries(curgap, curgap_num, gaps, gap_total, data)
        if (curgap_num != 1 and curgap[0] - gaps[curgap_num - 1][1] > 2) or (
            curgap_num == 1 and curgap[0] >= 2
        ):
            if (
                curgap_num != gap_total and gaps[curgap_num + 1][0] - curgap[1] > 2
            ) or (curgap_num == gap_total and len(data) - curgap[1] > 2):
                # 2-endpoint average matching
                m = (
                    (data[end[0]] + data[end[1]]) / 2
                    - (data[start[1]] + data[start[0]]) / 2
                ) / (xcoords[end[0]] - xcoords[start[1]])
                b = (data[start[1]] + data[start[0]]) / 2 - data[curgap[0]]
            elif curgap_num != gap_total or len(data) - curgap[1] == 2:
                # not 2 data points to match with on right side - endpoint matching w/ 2-pt average
                m = (data[end[0]] - (data[start[1]] + data[start[0]]) / 2) / (
                    xcoords[end[0]] - xcoords[start[1]]
                )
                b = (data[start[1]] + data[start[0]]) / 2 - data[curgap[0]]
            else:  # gap is at end of data
                m = 0
                b = (data[start[1]] + data[start[0]]) / 2 - data[curgap[0]]
        elif curgap_num != 1 or curgap[0] == 1:
            if (
                curgap_num != gap_total and gaps[curgap_num + 1][0] - curgap[1] > 2
            ) or (curgap_num == gap_total and len(data) - curgap[1] > 2):
                # not 2 data points to match with on left side - endpoint matching w/ 2-pt average
                m = ((data[end[0]] + data[end[1]]) / 2 - data[start[1]]) / (
                    xcoords[end[0]] - xcoords[start[1]]
                )
                b = data[start[1]] - data[curgap[0]]
            elif curgap_num != gap_total or len(data) - curgap[1] == 2:
                # not 2 data points to match with on either side - endpoint matching
                m = (data[end[0]] - data[start[1]]) / (
                    xcoords[end[0]] - xcoords[start[1]]
                )
                b = data[start[1]] - data[curgap[0]]
            else:  # gap is at end of data
                m = 0
                b = 0
        else:
            m = 0
            b = 0
    return m, b


def check_boundaries(curgap, curgap_num, gaps, gap_total, data):
    """Look for outliers at the current gap's boundaries.

    Parameters:
            curgap - tuple, start and end indices of the current gap
            curgap_num - int, index of current gap
            gaps - dict, stores all gaps in dataset
            gap_total - int, total number of gaps
            data - np.ndarray, the data being imputed
    """
    start_int = None
    start_ext = None
    end_int = None
    end_ext = None
    if (curgap_num != 1 and curgap[0] - gaps[curgap_num - 1][1] < 11) or curgap[0] < 11:
        start_ext = curgap[0] - 2
        start_int = curgap[0] - 1
    else:
        # find outliers using MAD (median average deviation) with 1st differences
        data_diff = []
        for i in range(curgap[0] - 11, curgap[0] - 1):
            data_diff.append(data[i + 1] - data[i])
        mad = median_abs_deviation(data_diff)

        start_ext = curgap[0] - 2
        start_int = curgap[0] - 1

        cutoff = 5 * mad
        outliers = [
            num > cutoff for num in data_diff
        ]  # positive difference indicates outlier at that index

        for i in range(0, 10):
            # outlier array in reverse order -> incrementing up outlier array decrements data index
            if outliers[i] or outliers[i + 1]:
                start_ext -= 1
                start_int -= 1
            else:
                break  # if no outliers, take indices and break for loop

    if (curgap_num != gap_total and gaps[curgap_num + 1][0] - curgap[1] < 11) or len(
        data
    ) - curgap[1] < 11:
        end_int = curgap[1] + 1
        end_ext = curgap[1] + 2
    else:
        # find outliers using MAD (median average deviation) with 1st differences
        data_diff = []
        for i in range(curgap[1] + 1, curgap[1] + 11):
            data_diff.append(data[i + 1] - data[i])
        mad = median_abs_deviation(data_diff)

        end_int = curgap[1] + 1
        end_ext = curgap[1] + 2

        cutoff = 5 * mad
        outliers = [
            num > cutoff for num in data_diff
        ]  # positive difference indicates outlier at that index

        for i in range(0, 10):
            # outlier array in reverse order -> incrementing up outlier array decrements data index
            if outliers[i] or outliers[i + 1]:
                end_int += 1
                end_ext += 1
            else:
                break  # if no outliers, take indices and break for loop

    if start_int < 0:
        start_int = None
    if start_ext < 0:
        start_ext = None
    if end_int > len(data) - 1:
        end_int = None
    if end_ext > len(data) - 1:
        end_ext = None

    start = (start_ext, start_int)
    end = (end_int, end_ext)
    return start, end


def fillgaps(datafile, method):
    """A method to fill gaps in data and preserve data noise characteristics.

    Parameters:
    datafile: a csv file formatted as follows: MJD, residual (us)
    method: a string indicating which of the 4 methods of data gap filling to use (reflect, reflect_invert, replica, replica_endpoint)
    """

    y = []
    x = np.array([])
    if datafile[-4:] == ".csv":  # if .csv file entered
        with open(datafile, "r") as f:
            lines = f.readline()
            while lines:
                vals = lines.split(",")
                x = np.append(x, float(vals[0]))
                y.append(float(vals[1]))
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

    fig, ax = plt.subplots(2, figsize=(12, 6))
    ax[0].errorbar(x, y, fmt=".")
    ax[0].grid(True)
    ax[0].set_xlabel("MJD")
    ax[0].set_ylabel("Residuals (us)")

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

    # impute gaps to the right of this run, then reverse and do left (4 methods)
    # i. Type 2 - reflect+invert
    if method == "reflect+invert":
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

    with open("result.csv", "w") as f:
        csvwriter = csv.writer(f)
        csvwriter.writerows(
            zip(xfilled, data)
        )  # write results into csv file in same format as input

    ax[1].errorbar(xfilled, data, fmt=".")
    ax[1].grid(True)
    ax[1].set_xlabel("MJD")
    ax[1].set_ylabel("Residuals (us)")

    plt.show()  # displays graphs of data before/after imputation


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError("Need two arguments: data file and imputation method.")

    fillgaps(sys.argv[1], sys.argv[2])
