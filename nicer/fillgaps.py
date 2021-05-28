"""
A script to fill data gaps in a way that retains the noise characteristic of the live data set and preserves AVAR curves. An implementation of the algorithm described in:
Howe and Schlossberger, "Strategy for Characterizing Frequency Stability Measurements having Multiple Dead Times"
"""

import argparse
import numpy as np
import copy


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


def fill(data, gaps, gap_total, gap_num, reverse=False):
    """A method to fill the gap_num-th gap in the data by reflecting and inverting data on one/both sides of the gap.

    Parameters:
            data - np.array, the data being imputed
            gaps - dict, stores the gap indeces with their respective number
            gap_total - int, total number of gaps
            gap_num - int, number of current gap
            reverse - bool, indicated whether dataset forward-oriented (gap being filled to the right of longest continuous data run) or reversed (left of longest run)
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
            data[curgap[0] : (curgap[0] + 1 + int(size / 2))] = left_pts_to_reflect[
                ::-1
            ]
            data[(curgap[1] - int(size / 2)) : curgap[1] + 1] = right_pts_to_reflect[
                ::-1
            ]
            # calculate linear slope to apply to fill pts
            m, b = calculate_slope(
                data, curgap, gap_num, gaps, gap_total, len(pts_to_reflect)
            )
            for j in range(0, size):
                data[curgap[0] + j] = m * j + data[curgap[0] + j] + b
            # if odd number of spaces in gap, average first left and last right values and put in middle gap value
            if (size / 2) % 2 == 0.5:
                mid_val = (left_pts_to_reflect[0] + right_pts_to_reflect[-1]) / 2
                data[curgap[1] - int(size / 2)] = mid_val
            if gap_num == gap_total:
                gaps.pop(gap_num)
                gap_num = gap_num - 1
            # remove gap from gap map
            for k in range(gap_num, gap_total):
                gaps[k] = gaps.pop(k + 1)
            gap_total = gap_total - 1
        else:
            # invert and add linear slope to match gap-fill ends to existing data
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
            # calculate linear slope to apply to fill pts
            data[curgap[0] : (curgap[1] + 1)] = pts_to_reflect[::-1]
            m, b = calculate_slope(
                data, curgap, gap_num, gaps, gap_total, len(pts_to_reflect)
            )
            for j in range(0, size):
                data[curgap[0] + j] = m * j + data[curgap[0] + j] + b
            if gap_num == gap_total:
                gaps.pop(gap_num)
                gap_num = gap_num - 1
            else:
                for k in range(gap_num, gap_total):
                    gaps[k] = gaps.pop(k + 1)
            gap_total = gap_total - 1
    return gap_num, gap_total


def calculate_slope(data, curgap, curgap_num, gaps, gap_total, fill_pts_len):
    """A method to calculate a linear slope to add to the fill points for the gap and ensure 2 end-point matching when possible.

    Parameters:
            data - np.array, the data being imputed
            curgap - tuple, (start, end) indeces for current gap
            curgap_num - int, number of current gap
            gaps - dict, stores the gap indeces with their respective number
            gap_total - int, total number of gaps
            fill_pts_len - int, the length of the data to fill the gap
    """
    m = 0
    b = 0
    if fill_pts_len > 1:
        if (curgap_num != 1 and curgap[0] - gaps[curgap_num - 1][1] > 2) or (
            curgap_num == 1 and curgap[0] >= 2
        ):
            if (
                curgap_num != gap_total and gaps[curgap_num + 1][0] - curgap[1] > 2
            ) or (curgap_num == gap_total and len(data) - curgap[1] > 2):
                # 2-endpoint average matching
                m = (
                    (data[curgap[1] + 1] + data[curgap[1] + 2]) / 2
                    - (data[curgap[0] - 1] + data[curgap[0] - 2]) / 2
                ) / (fill_pts_len + 4)
                b = (data[curgap[0] - 1] + data[curgap[0] - 2]) / 2 - data[curgap[0]]
            elif curgap_num != gap_total or len(data) - curgap[1] == 2:
                # not 2 data points to match with on right side - endpoint matching w/ 2-pt average
                m = (
                    data[curgap[1] + 1]
                    - (data[curgap[0] - 1] + data[curgap[0] - 2]) / 2
                ) / (fill_pts_len + 3)
                b = (data[curgap[0] - 1] + data[curgap[0] - 2]) / 2 - data[curgap[0]]
            else:  # gap is at end of data
                m = 0
                b = (data[curgap[0] - 1] + data[curgap[0] - 2]) / 2 - data[curgap[0]]
        elif curgap_num != 1 or curgap[0] == 1:
            if (
                curgap_num != gap_total and gaps[curgap_num + 1][0] - curgap[1] > 2
            ) or (curgap_num == gap_total and len(data) - curgap[1] > 2):
                # not 2 data points to match with on left side - endpoint matching w/ 2-pt average
                m = (
                    (data[curgap[1] + 1] + data[curgap[1] + 2]) / 2
                    - data[curgap[0] - 1]
                ) / (fill_pts_len + 2)
                b = data[curgap[0] - 1] - data[curgap[0]]
            elif curgap_num != gap_total or len(data) - curgap[1] == 2:
                # not 2 data points to match with on either side - endpoint matching
                m = (data[curgap[1] + 1] - data[curgap[0] - 1]) / (fill_pts_len + 2)
                b = data[curgap[0] - 1] - data[curgap[0]]
            else:  # gap is at end of data
                m = 0
                b = 0
        else:
            m = 0
            b = 0
    return m, b


def fillgaps(data, method):
    """A method to fill gaps in data and preserve data noise characteristics.

    Parameters:
    data: an array of data with gaps represented as numpy.nan
    method: a string indicating which of the 4 methods of data gap filling to use (reflect, reflect_invert, replica, replica_endpoint)
    """

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
                gap_to_impute, gap_total = fill(data, gaps, gap_total, gap_to_impute)
                if gap_to_impute is None:
                    return

            # spacially reverse dataset, continue until beginning reached
            for i in range(gap_to_impute, 0, -1):
                gap_to_impute, gap_total = fill(data, gaps, gap_total, i, reverse=True)
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
