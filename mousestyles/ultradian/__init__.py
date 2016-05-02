from __future__ import print_function, absolute_import, division

import numpy as np
import pandas as pd

import mousestyles.data as data


def aggegate_interval(strain, mouse, feature, bin_width):
    """
    data loaded from data.load_intervals(feature)

    Parameters
    ---------------
    feature: {"AS", "Food", "IS", "M_AS", "M_IS", "Water",
              "Distance", "AS_Intensity", "AS_prob"}
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
                time_behaviour[
                    start_index + j *
                    bin_count] += bin_length * end_index - start_time
                time_behaviour[end_index + j *
                               bin_count] += end_time % bin_length
            else:
                time_behaviour[
                    start_index + j *
                    bin_count] += bin_length * (start_index + 1) - start_time
                time_behaviour[end_index + j *
                               bin_count] += end_time % bin_length
                time_behaviour[start_index + j * bin_count +
                               1:end_index + j * bin_count] += bin_length

    if feature == 'F' or feature == 'W':
        all_feature = data.load_all_features()
        group = all_feature[
            ["strain", "mouse", "day", "hour", "Food", "Water"]].groupby(
            ["strain", "mouse", "day"]).sum()
        group = group.reset_index()
        group['day'] -= 5
        mouse_data = group.loc[(group['strain'] == strain) &
                               (group['mouse'] == mouse)]
        for i in mouse_data['day'].astype('int'):
            if feature == 'F':
                food_amount = float(mouse_data['Food'][mouse_data['day'] == i])
                time_behaviour[
                    (bin_count * i):(bin_count * (i + 1) - 1)] /= sum(
                    time_behaviour[(bin_count * i):(bin_count * (i + 1) - 1)])
                time_behaviour[(bin_count * i):(bin_count *
                                                (i + 1) - 1)] *= food_amount
            else:
                food_amount = float(mouse_data['Water'][
                                    mouse_data['day'] == i])
                time_behaviour[
                    (bin_count * i):(bin_count * (i + 1) - 1)] /= sum(
                    time_behaviour[(bin_count * i):(bin_count * (i + 1) - 1)])
                time_behaviour[(bin_count * i):(bin_count *
                                                (i + 1) - 1)] *= food_amount

    ts = pd.Series(time_behaviour, index=pd.date_range(
        '01/01/2014', periods=len(time_behaviour),
        freq=str(bin_width) + 'min'))
    return(ts)


def aggregate_movement(strain, mouse, bin_width):
    """
    data loaded from data.load_movement(strain, mouse, days = ...)

    Parameters
    ---------------
    bin_width: number of minutes of time interval for data aggregation

    Returns
    ----------
    ts: pandas.tseries
        a pandas time series of length (#day)*24(hour)*60(minute)/n
    """
    # determine number of days
    intervals = data.load_intervals('IS')
    mouse_data = intervals.loc[
        (intervals['strain'] == strain) & (intervals['mouse'] == mouse)]
    days = sorted(np.unique(mouse_data['day']))

    # build data frame
    bin_count = int(24 * 60 / bin_width)
    time_movements = np.repeat(0.0, bin_count * len(days))
    bin_length = bin_width * 60
    for j in days:
        M = data.load_movement(strain, mouse, day=int(j))
        distance_df = pd.DataFrame({"start": M["t"].values[0:-1], "end": M["t"].values[1:],
                                    "distance": np.linalg.norm(M[["x", "y"]].values[1:] -
                                                               M[["x", "y"]].values[0:-1], axis=1)})
        start_end = data.load_start_time_end_time(strain, mouse, j)
        start = np.asarray(distance_df['start']) - start_end[0]
        end = np.asarray(distance_df['end']) - start_end[0]
        dist = distance_df['distance']
        for i in range(len(start)):
            start_time = start[i]
            end_time = end[i]
            start_index = int(start_time / (bin_width * 60))
            end_index = int(end_time / (bin_width * 60))
            if start_index == end_index:
                time_movements[start_index + j *
                               bin_count] += dist[i]
            else:
                time_movements[end_index + j *
                               bin_count] += end_time % bin_length / (end_time - start_time) * dist[i]
                time_movements[
                    start_index + j * bin_count] += dist[i] - end_time % bin_length / (end_time - start_time) * dist[i]
    ts = pd.Series(time_movements, index=pd.date_range(
        '01/01/2014', periods=len(time_movements),
        freq=str(bin_width) + 'min'))
    return(ts)
