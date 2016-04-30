from __future__ import print_function, absolute_import, division

import numpy as np
import pandas as pd

import datetime
import os

import mousestyles.data as data


def aggegate_interval(strain, mouse, feature, bin_width):
    """
    data loaded from data.load_intervals(feature)

    Parameters
    ---------------
    feature: {"AS", "Food", "IS", "M_AS", "M_IS", "Water", "Distance", "AS_Intensity", "AS_prob"}
    bin_width: number of minutes of time interval for data aggregation

    Returns
    ----------
    ts: pandas.tseries
        a pandas time series of length 12(day)*24(hour)*60(minute)/n
    """
    # load data
    intervals = data.load_intervals(feature)
    mouse_data = intervals.loc[
        (intervals['strain'] == strain) & (intervals['mouse'] == mouse)]

    # build data frame
    days = sorted(np.unique(mouse_data['day']))
    bin_count = int(24 * 60 / bin_width)
    time_behaviour = np.repeat(0.0, bin_count * len(days))
    bin_length = bin_width * 60

    for j in days:
        df = mouse_data.loc[mouse_data['day'] == j]
        start_end = data.load_start_time_end_time(strain, mouse, j)
        start = np.asarray(df['start']) - start_end[0]
        end = np.asarray(df['stop']) - start_end[0]

        for i in range(len(start)):
            start_time = start[i]
            end_time = end[i]
            start_index = int(start_time / (bin_width * 60))
            end_index = int(end_time / (bin_width * 60))
            if start_index == end_index:
                time_behaviour[start_index + j *
                               bin_count] += end_time - start_time
            elif end_index - start_index == 1:
                time_behaviour[start_index + j *
                               bin_count] += bin_length * end_index - start_time
                time_behaviour[end_index + j *
                               bin_count] += end_time % bin_length
            else:
                time_behaviour[
                    start_index + j * bin_count] += bin_length * (start_index + 1) - start_time
                time_behaviour[end_index + j *
                               bin_count] += end_time % bin_length
                time_behaviour[start_index + j * bin_count +
                               1:end_index + j * bin_count] += bin_length
    ts = pd.Series(time_behaviour, index=pd.date_range(
        '01/01/2014', periods=len(time_behaviour), freq=str(bin_width) + 'min'))
    return(ts)
