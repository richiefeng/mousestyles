from __future__ import print_function, absolute_import, division

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

import mousestyles.data as data


def aggegate_interval(strain, mouse, feature, bin_width):
    """
    Aggregate the interval data based on n-minute time
    intervals, return a time series.

    Parameters
    ----------

    strain: int
        nonnegative integer indicating the strain number
    mouse: int
        nonnegative integer indicating the mouse number
    feature: {"AS", "F", "IS", "M_AS", "M_IS", "W", "AS_prob"}
    bin_width: number of minutes of time interval for data aggregation

    Returns
    -------

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
    Aggregate the movement data based on n-minute
    time intervals, return a time series.

    Parameters
    ----------
    strain: int
        nonnegative integer indicating the strain number
    mouse: int
        nonnegative integer indicating the mouse number
    bin_width: number of minutes of time interval for data aggregation

    Returns
    -------
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
        distance_df = pd.DataFrame({"start": M["t"].values[0:-1],
                                    "end": M["t"].values[1:],
                                    "distance":
                                    np.linalg.norm(M[["x", "y"]].values[1:] -
                                                   M[["x", "y"]].values[0:-1],
                                                   axis=1)})
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
                time_movements[
                    end_index + j * bin_count] += end_time % \
                    bin_length / (end_time - start_time) * dist[i]
                time_movements[
                    start_index + j * bin_count] += dist[i] - \
                    end_time % bin_length / (end_time - start_time) * dist[i]
    ts = pd.Series(time_movements, index=pd.date_range(
        '01/01/2014', periods=len(time_movements),
        freq=str(bin_width) + 'min'))
    return(ts)


def seasonal_decomposition(strain, mouse, feature, bin_width, period_length):
    """
    Apply seasonal decomposition model on the time series
    of specified strain, mouse, feature and bin_width.

    Parameters
    ----------
    strain: int
        nonnegative integer indicating the strain number
    mouse: int
        nonnegative integer indicating the mouse number
    feature: {"AS", "F", "IS", "M_AS", "M_IS", "W"}
    bin_width: int
        number of minutes, the time interval for data aggregation
    period_length: int
        number of hours, usually the significant period
        length indicated by Lomb-scargle model

    Returns
    -------
    res: statsmodel seasonal decomposition object
        seasonal decomposition result for the mouse.
        Check the seasonal decomposition plot by res.plot(),
        seasonl term and trend term by res.seasonal and
        res.trend separately.

    Examples
    --------
    >>> res = seasonal_decomposition(strain=0, mouse=0, feature="W",
                                     bin_width=30, period_length = 24)
    """
    freq = int(period_length * 60 / bin_width)
    ts = aggegate_interval(strain=strain, mouse=mouse,
                           feature=feature, bin_width=bin_width)
    res = sm.tsa.seasonal_decompose(ts.values, freq=freq, model="additive")
    return(res)


def strain_seasonal(strain, mouse, feature, bin_width, period_length):
    """
    Use seansonal decomposition model on the time series
    of specified strain, mouse, feature and bin_width.
    return the seasonal term and the plot of seasonal term
    by mouse of a set of mouses in a strain

    Parameters
    ----------
    strain: int
        nonnegative integer indicating the strain number
    mouse: list, set or tuple
        nonnegative integer indicating the mouse number
    feature: {"AS", "F", "IS", "M_AS", "M_IS", "W"}
    bin_width: int
        number of minutes, the time interval for data aggregation
    period_length: int
        number of hours, usually the significant period
        length indicated by Lomb-scargle model

    Returns
    -------
    seasonal_all: numpy array containing the seasonal term for every
        mouse indicated by the input parameter
    seasonal_plot: plot of seasonal term by mouse

    Examples
    --------
    >>> res = strain_seasonal(strain=0, mouse={0, 1, 2, 3}, feature="W",
                              bin_width=30, period_length = 24)
    """
    # seasonal decomposition
    seasonal_all = np.array([])
    freq = int(period_length * 60 / bin_width)
    for m in mouse:
        res = seasonal_decomposition(
            strain, m, feature, bin_width, period_length)
        seasonal_all = np.append(seasonal_all, res.seasonal[0:freq])
    seasonal_all = seasonal_all.reshape([len(mouse), -1])
    # seasonal plot
    seasonal_plot = plt.figure()
    ax = seasonal_plot.add_subplot(1, 1, 1)
    time = np.arange(0, period_length, bin_width / 60)
    for i in np.arange(len(mouse)):
        ax.plot(time, seasonal_all[i, :])
    ax.legend(['mouse' + str(i) for i in np.arange(len(mouse))],
              loc='upper right', prop={'size': 8})
    plt.xlabel('Time')
    plt.xlabel('Time')
    plt.ylabel('Seasonal Term')
    return seasonal_all, seasonal_plot
