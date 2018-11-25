# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 12:38:41 2018

@author: JEROMLU2
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.gridspec as gridspec

import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures, minmax_scale
from sklearn.model_selection import cross_val_score, LeaveOneOut, StratifiedKFold
from sklearn.metrics import make_scorer, r2_score


def create_RS_features(X_df):
    '''
    It creates interactions and quadratic terms - defined by degree = 2.
    
    Parameters
    ----------
    X_df : pd.DataFrame
        data frame of only main features
    
    Returns
    -------
    X_df_RS (pd.DataFrame)
    
    '''
    
    poly = PolynomialFeatures(2)
    poly.fit(X_df)
    
    #create dataframe properly named, of full response surface (second order)
    X_df_RS = pd.DataFrame(poly.transform(X_df),
                                 columns = poly.get_feature_names(X_df.columns),
                                 index = X_df.index)
    return X_df_RS
    
def press_statistic(y_true, y_pred, xs):
    """
    Calculation of the `Press Statistics <https://www.otexts.org/1580>`_
    
    Parameters
    ----------
    y_true : pd.Series
        pandas series of measured targets
    y_pred : pd.Series
        pandas series of predicted values
    xs : pd.DataFrame
    
    Returns
    -------
    calculated PRESS value (float)
    """
    res = y_pred - y_true
    hat = xs.dot(np.linalg.pinv(xs))
    den = (1 - np.diagonal(hat))
    sqr = np.square(res/den)
    return sqr.sum()

def PRESS(y_true, y_pred, X):
    
    res = y_pred - y_true
    hat = X.dot(np.linalg.inv(X.T.dot(X)).dot(X.T))
    den = (1 - np.diagonal(hat))
    sqr = np.square(res/den)
    return sqr.sum()
    

def predicted_r2(y_true, y_pred, xs = None):
    """
    Calculation of the `Predicted R-squared <https://rpubs.com/RatherBit/102428>`_
    TODO: validate this statement: as in MODDE software
    
        Parameters
    ----------
    y_true : pd.Series
        pandas series of measured targets
    y_pred : pd.Series
        pandas series of predicted values
    xs : pd.DataFrame
        data frame of feature values
    
    Returns
    -------
    calculated PRESS value (float)
    """
    press = press_statistic(y_true=y_true,
                            y_pred=y_pred,
                            xs=xs
    )
    
    press2 = PRESS(y_true=y_true,
                            y_pred=y_pred,
                            X=xs
    )

    sst  = np.square( y_true - y_true.mean() ).sum()
    
    return 1 - press / sst
 
def r2(y_true, y_pred):
    """
    Calculation of the unadjusted r-squared, goodness of fit metric
    tega ne rabom, ker je ze vgrajeno v scikit learn package
    vseeno ohranim, za formulo
    """
    sse  = np.square( y_pred - y_true ).sum()
    sst  = np.square( y_true - y_true.mean() ).sum()
    return 1 - sse/sst

def stepwise_feat_selection(y_df, X_df, tresh = 0.05, method = 'p-value_tresh',
                            cv = 5, restrict = False, main_effects = None,
                            add_steps = 4):
    '''
    Stepwise features selection (elimination). We start with full response surface
    model and in each step the feature with maximum p-value is eliminated.
    
    Parameters
    ----------
    y_df : pd.Series
        vectore of targets (responses)
    X_df : pd.DataFrame
        matrix of features (factors), second dimension (number of columns) 
        has to be same as length of y_df, The feature matrix has to include 
        constant term.
    tresh : float 
        floating point number in the interval (0 and 1). This is used as 
        threshold for feature elimination.
    method : string
        TODO
    cv : int
        number of subgroups. Should not be larger than the number of samples
        (rows in X_df)
    restrict : bool
        if True main features are not removed if they are present in 
        any of other combinations
    add_steps : int
        for how many steps to continue, it could happen that you run out
        of the features WATCH it!!!
        
    
    Returns
    -------
    results (statsmodels), model (statsmodels), 
    history (list of results), Q_scores (list of cross valiadtion scores)
    
    '''
    #some initial settings
    add_steps = add_steps
    max_iterations = X_df.shape[-1] - 2
    
    #counter for number of loops
    i = 0
    i_final = max_iterations + add_steps
    if main_effects is None:
        main_effects = [effect for effect in X_df.columns if ('^' not in effect) and (' ' not in effect)]
    #print(main_effects,'\n\n')
    
    #initialization of the model
    model = sm.OLS(y_df, X_df)
    results = model.fit()
    
    
    model = LinearRegression(fit_intercept = False) #since I added constant term by hands
    
    Q2_score = make_scorer(predicted_r2, greater_is_better=True, xs = X_df.values)
    Q2_res = Q2_score(LinearRegression().fit(X_df, y_df), X_df, y_df)
    scores = cross_val_score(model, X_df, y_df, cv = cv, scoring = 'r2') #r2

    
    #initialization of history and Q_scores
    history = [results] #statsmodels
    cross_scores = [scores] #scikit-learn
    Q2_scores = [Q2_res]
    
    p_values = results.pvalues
    p_values_temp = p_values.drop('1').copy()
    if restrict == True:
        rest_of_mains = p_values.drop(main_effects, errors = 'ignore')
        #print(rest_of_mains, '\n')
        for name in rest_of_mains.index:
            #print('name ', name)
            for main_effect in main_effects:
                #print('effect' ,main_effect)
                if (main_effect in name) and (main_effect in p_values_temp.index):
                    p_values_temp.drop(main_effect, inplace= True)


    droped_name = p_values_temp.idxmax()
    names_to_keep = p_values.drop(droped_name).index
    
    X_df_new = X_df[names_to_keep]
    
    #loop_stop_condition could as well be: 
    #.......p_values.max() > tresh
    #with this condition actually last model is not the best
    print(max_iterations)
    while (i - add_steps < i_final) and (i < max_iterations):
        if len(p_values) == 1:
            break
        #fit model with the new data
        #statsmodels
        model = sm.OLS(y_df, X_df_new)
        results = model.fit()
        p_values = results.pvalues
        p_values_temp = p_values.drop('1').copy()
        

        
        
        #sickit learn for cross validation
        model = LinearRegression(fit_intercept = False) #since I added constant term by hands
        scores = cross_val_score(model, X_df_new, y_df, cv=cv, scoring = 'r2')
        Q2 = Q2_score(LinearRegression().fit(X_df_new, y_df), X_df_new, y_df)
        
        #creation of history
        history.append(results) #statsmodels
        cross_scores.append(scores) #scikit-learn
        Q2_scores.append(Q2) #my function + sklearn
        
        if restrict == True:
            rest_of_mains = p_values.drop(main_effects, errors = 'ignore')
            #print(rest_of_mains, '\n')
            for name in rest_of_mains.index:
                #print('name ', name)
                for main_effect in main_effects:
                    #print('effect' ,main_effect)
                    if (main_effect in name) and (main_effect in p_values_temp.index):
                        p_values_temp.drop(main_effect, inplace= True)
        #print(p_values_temp)
        droped_name = p_values_temp.idxmax()
        names_to_keep = p_values.drop(droped_name).index
        #print(names_to_keep, '\n\n')
        X_df_new = X_df_new[names_to_keep]
        
        print(i)
        #number of loops counter
        i = i + 1
        
        #stop condition
        if (p_values_temp.max() < tresh) and i_final > max_iterations:
            i_final = i
        
        #reporting progress
        text = 'Step: {0}: {1:.<25} was removed with p value of {2:5.3g}'
        print(text.format(i, droped_name, max(p_values_temp)))
    print('Number of steps: ', i)
    print('Best model at step: ', i_final)
    
    
    #final fit of the model
    X_df_final = X_df[history[i_final].params.index]
    model = sm.OLS(y_df, X_df_final)
    results = model.fit()
    
    #show the final results
    print(results.summary())
    
    return results, model, history, cross_scores, Q2_scores


def code_df(df, df_upper, df_lower):
    '''
    Codes the data frame object so that there are only ones and zeros present in it.
    Parameters
    -----------
    df ... the feature matrix(pandas data frame object) that should be transformed
    df_upper ... pandas series with upper values
    df_lower ... pandas series with lower values
    
    Returns
    --------
    coded data frame (pd.DataFrame)
    '''
    cols = df.columns
    #same as below
    pd.DataFrame(minmax_scale(df, (-1, 1)), columns = df.columns, index = df.index)
    if set(cols).issubset(df_upper.index) and set(cols).issubset(df_lower.index):
        return (df -((df_upper[cols]+df_lower[cols])/2))/(df_upper[cols]-df_lower[cols])*2
    else:
        print('Columns of df and upper and lower do not agree')
        return pd.DataFrame([])
    
def scatter_features_target_pot(X_df, y_df, target_name = 'target'):
    '''
    Plots targets in dependance for each fature that is present in X_df
    Parameters
    ----------
    X_df : pd.DataFrame
        data frame of all feature settings, should contain only main fatures
    y_df : pd.Series
        series of targets
    target_naem : str
        name of target written to y axis
    '''
    for feature in X_df.columns:
        fig, ax = plt.subplots(figsize = (7, 4))
        ax.scatter(X_df[feature], y_df)
        ax.set_xlabel(feature)
        ax.set_ylabel(target_name)

def scatter_plot_QA(y, features, mark_center = True, center_points = None):
    '''
    Scater plot of all results (responses). It marks the center point runs with different color
    '''
    
    X = features
    if center_points is None:       
        center_points_mask = X.eq(0.0, axis=0).all(axis=1)
        center_points = X[center_points_mask].index
        print('Center points are:\n',list(center_points))
    colors = pd.Series('k', index = X.index)
    colors.loc[center_points] = 'r'
    fig, ax = plt.subplots(figsize = (10, 4))
    ax.scatter(np.arange(len(y)), y, c = colors)
    ax.set_title('Scatter plot of targets data points ')
    ax.set_ylabel(y.name)
    fig.tight_layout()
    plt.show()
    return ax

def plot_coeff_effects(results, sorted_by = 'pvalues'):
    '''
    Plots horizontal bars of logarithem of p values and effects.
    Bars in plots are sorted by their p values.
    '''
    bar_width = 0.5
    fig, ax = plt.subplots(1,2, figsize = (11, 5))
    res = pd.DataFrame([])
    res['params'] = results.params
    res['pvalues'] = results.pvalues
    const = res.iloc[:1]
    rest = abs(res.iloc[1:]).sort_values(sorted_by, ascending  = False)
    order = list(rest.index) + ['const']
    ax[0].barh(rest.index, rest.params, height = bar_width)
    ax[0].barh(const.index, abs(const.params), height = bar_width)
    ax[0].set_yticklabels(order)
    ax[0].set_title('Effects')
    order = list(rest.index) + ['1']
    ax[1].barh(rest.index, -np.log(rest.pvalues), height = bar_width)
    ax[1].barh(const.index, -np.log(const.pvalues), height = bar_width)
    ax[1].set_yticklabels('')
    ax[1].set_title('Log worth of effects.')
    
    trans = blended_transform_factory(ax[1].transData, ax[1].transAxes)
    ax[1].vlines(-np.log(0.05),0 , 1, transform = trans)
    fig.tight_layout()
    plt.show()
    return ax

def plot_history_of_stepwise(target_name, history, cross_scores, Q2_scores):
    '''
    Plots history of the fit (AIC and BIC for each result).
    
    Parameeters
    -----------
    target_name : str
        name of the target (response) that is under investigation
    history : list
        list of statsmodels results
    cross_scores : list
        list of all cross validation scores
    Q2_scores : list
        list of all Q2 scores calculated from PRESS - see MODDE user guide
        
    Returns
    -------
    ax (mpl.axes)
    '''
    
    all_factors = len(history[0].params)
    history_df = pd.Series(history)
    fig, ax_list = plt.subplots(1, 2, figsize = (14, 5))
    ax = ax_list[0]
    func = lambda x : x.bic
    y_vls = history_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls,y_vls,'-o',label = 'bic')
    func = lambda x : x.aic
    y_vls = history_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls,y_vls,'-o', label =  'aic')
    ax.set_xticks(x_vls)
    x_tick_labels = [ str(el) + '\n' + str(i) for i, el in enumerate(x_vls)]
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel('Number of factors')
    ax.set_title(target_name)
    ax.legend()
    
    history_df = pd.Series(history)
    cross_scores_df = pd.Series(cross_scores)
    
    ax = ax_list[1]
    func = lambda x : x.rsquared
    y_vls = history_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls ,y_vls,'-o',label = 'R^2')
    func = lambda x : x.rsquared_adj
    y_vls = history_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls,y_vls,'-o', label =  'R^2_adj')
    func = lambda x : np.mean(x)
    y_vls = cross_scores_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls,y_vls,'-o', label =  'cross_r2_mean')
    func = lambda x : np.max(x)
    y_vls = cross_scores_df.apply(func)
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls,y_vls,'-o', label =  'cross_r2_max')
    y_vls = Q2_scores
    x_vls = all_factors - np.arange(len(y_vls))
    ax.plot(x_vls, y_vls,'-o', label =  'predicted_R2')
    
    ax.set_xlabel('Number of factors')
    ax.set_xticks(x_vls)
    x_tick_labels = [ str(el) + '\n' + str(i) for i, el in enumerate(x_vls)]
    ax.set_xticklabels(x_tick_labels)
    ax.set_ylim(bottom=0)
    ax.set_title(target_name)
    ax.legend()
    fig.tight_layout()
    plt.show()
    return ax



def recalculate_feature_matrix(results, feat_settings, x_var, y_var, X, Y, columns = None):  
    if columns is not None:
        exog = pd.DataFrame(columns = columns)
    else:
        exog = pd.DataFrame()  
    exog[x_var] = X.ravel()
    exog[y_var] = Y.ravel()
    for key in feat_settings:
        exog[key] = feat_settings[key]

    exog = create_RS_features(exog)
    return exog[results.params.index]

def recalculate_1D(results, x_var, x, feat_settings, columns):
    '''
    Recalculates the feature matrix for the prediction. 
    Only one column changes all others are kept constant (1D).
    
    Parameters
    ----------
    results : statsmodels result object
    x_var : str
    x : np.array
    feat_settings : dict
    columns :  
        
    Returns
    -------
    exog (pd.DataFrame)
        
    '''

    if columns is not None:
        exog = pd.DataFrame(columns = columns)
    else:
        exog = pd.DataFrame()
    exog[x_var] = x
    for feat in feat_settings:
        exog[feat] = feat_settings[feat]
    
    exog = create_RS_features(exog)
    return exog[results.params.index]


def get_main_effects(tab, const = False):
    '''
    Function that tries to determine which are main effects.
    There should be no white spaces in the names of main effects, otherwise this will not work properly.
    
    '''
    main_effects = [effect for effect in tab if ('^' not in effect) and (' ' not in effect)]
    if const:
        return main_effects[1:]
    else:
        return main_effects

def plot_2D_targets(ax, data):
    pass
    
def plot_1D_target(ax, x_var, data):
    '''
    Plot one profile.
    '''
    x, y = data
    curve = ax.plot(x, y, label = x_var)
    ax.set_title(x_var)
    return curve

def prediction_profiler(results, feat_limits, columns):
    '''
    Plots interactive prediction profiler. Not interactive yet.
    TODO: interactive
    
    Parameters
    ----------
    results : statsmodels result object
        results of the statsmodels fit (RegressioReslt object). I use .predict()
        method of this object to predict dependent variable
        as well as some other information that is kept in it. I also extract 
        the dependent variable name from the result object.
    feat_limits : dict
        dictonary where the limits of the plot are set for each of the independent 
        variables. It should contain all the independent variables.
    colums : list
        list of all the main features that were used before the start of the stepwise
        elimination of the features
        
    Returns
    -------
    fig (mpl.figure), ax_list (list of mpl.axes), sldr (list of mpl sliders)
    '''
    #get target (response) name from the results object, this is quite specific
    target_name = _get_target_name(results)    
    #get all possible main effects - no interactions and qudratic terms
    main_effects = get_main_effects(columns)
    
    #number of points to plot
    n_points = 20
    
    #get onls relevant effects that are present in result object
    rel_effects = get_main_effects(results.params.index, True)
    n_feat = len(rel_effects)
    
    #set scale and with of the sample according to the number of relevant features
    scale = 0.8
    wid = n_feat * 3.4 * scale
    
    
    fig = plt.figure(figsize= (wid, 4 * scale))
    ax_list = []
    
    gs = gridspec.GridSpec(2, n_feat,
                           height_ratios=[10, 1]
                           )
    
    for g in gs:
            ax_list.append(plt.subplot(g))
    ax_list[0].set_ylabel(target_name)
    ax_sldr = {}
    sldr = {}
    lns = {}
    axcolor = 'lightgoldenrodyellow'

    y_max = []
    y_min = []
    #for each of the relevant main effects we calculate and plot 1D graph
    for i, feature in enumerate(rel_effects):

        #set the independant variable name
        x_var = rel_effects[i]
        
        #reset the main effects
        main_effects = get_main_effects(columns)
        main_effects.remove(x_var)
        
        #set the features settings - this should be variable not all zero
        feat_settings = {}
        for feature in main_effects:
            feat_settings[feature] = 0
        
        #calculate the predictions and save it in z variable
        x = np.linspace(-1, 1, n_points)
        exog = recalculate_1D(results, x_var, x, feat_settings, columns)
        z = results.predict(exog)
        y_min.append(z.min())
        y_max.append(z.max())

        #rescale the data
        lower, upper = feat_limits[x_var]
        x = (upper + lower) / 2 +  x * (upper - lower) /2
    
        ln, = plot_1D_target(ax_list[i], x_var, [x, z])
        lns[x_var] = ln
        ax_sldr[x_var] = ax_list[i + n_feat]
        ax_sldr[x_var].set_facecolor(axcolor)
        #valinit = feat_settings[x_var]
        sldr[x_var] = Slider(ax_sldr[x_var], '', -1, 1, valinit = 0)
        #sldr[x_var].on_changed(update)
    
    
    def update(value):
        y_max = []
        y_min = [] 
        for i, feature in enumerate(rel_effects):
           
            #set the independant variable name
            x_var = rel_effects[i]
            
            #reset the main effects
            main_effects = get_main_effects(columns)
            main_effects.remove(x_var)
            
            #set the features settings according to the sliders
            feat_settings = {}
            for feature in main_effects:
                if feature in sldr:
                    feat_settings[feature] = sldr[feature].val
                else:
                    feat_settings[feature] = 0
            
            #calculate the predictions and save it in z variable
            x = np.linspace(-1, 1, n_points)
            exog = recalculate_1D(results, x_var, x, feat_settings, columns)
            z = results.predict(exog)
            y_min.append(z.min())
            y_max.append(z.max())
            lns[x_var].set_ydata(z)
        
        for i, feature in enumerate(rel_effects):
            ax_list[i].set_ylim(min(y_min), max(y_max))
        fig.canvas.draw_idle()


    #connect the sliders with the change effect
    for key in sldr:
        sldr[key].on_changed(update)
    for i, feature in enumerate(rel_effects):
        ax_list[i].set_ylim(min(y_min), max(y_max))    
    
    fig.tight_layout()
    #plt.show()
    #have to return sldr otherwise ti looses the reference and sliders are not responsive!!!
    return fig, ax_list, sldr

def _get_target_name(results):
    '''
    It gets the target (response) name from the statsmodel RegressionResults object
    it could be version specific
    '''
    return results.summary().tables[0].data[0][1]
        
def prediction_contour_plot(results, x_var, y_var, feat_limits, c_extent, columns):
    '''
    Plots 2D filled contour plot of one target (response) in dependence on two 
    features x_var and y_var. For rest of the main effects in the result object 
    it creates sliders with which one may control the settings of other features.
    
    Parameters
    ----------
    results : statsmodels result object
        results of the statsmodels fit (RegressioReslt object). I use .predict()
        method of this object to predict dependent variable
        as well as some other information that is kept in it. I also extract 
        the dependent variable name from the result object.
    x_var : str
        name of the independent variable (feature, factor) 
        which will be plotted on x axis
    y_var : str
        name of the independent variable (feature, factor) 
        which will be plotted on y axis
    feat_limits : dict or pd.DataFrame
        dictonary where the limits of the plot are set for each of the independent 
        variables. It should contain all the independent variables. Two columns 
        'lower ' and 'upper' that contain upper and lower limits of the data.
    c_extent : list
    colums : list
        list of all the main features that were used before the start of the stepwise
        elimination of the features.
        
    Returns
    -------
    Nothing
    
    '''
    #get target (response) name from the results object, this is quite specific
    target_name = _get_target_name(results)
    
    #get main effects from the columns, required for
    main_effects = get_main_effects(columns)
    #remove tha names which are ploted
    main_effects.remove(x_var)
    main_effects.remove(y_var)
    
    #number of points to plot, total size n_points^2
    n_points = 10
    #c_extent = code_df(extent, X_upper_lim, X_lower_lim).values.T.ravel()
    #creating number of points
    x = np.linspace(c_extent[0], c_extent[1], n_points)
    y = np.linspace(c_extent[2], c_extent[3], n_points)
    X, Y = np.meshgrid(x, y) # grid of points

    #creation of figure with proper size and remove the grid
    fig, ax_im = plt.subplots(figsize = (10, 9))
    ax_im.grid(False)
    
    #prepare data
    feat_settings = {}
    for feature in main_effects:
        feat_settings[feature] = 0
    #reacalculate feature matrix
    exog = recalculate_feature_matrix(results, feat_settings, x_var, y_var, X, Y, columns);
    
    #predict the values and write data to Z variable
    predicted = results.predict(exog)
    Z = predicted.values.reshape(X.shape)
    
    #calculate true extent
    extent = pd.Series([])
    
    n = (feat_limits.loc['lower'] + feat_limits.loc['upper']) / 2
    k = (feat_limits.loc['upper'] - feat_limits.loc['lower']) / 2
    
 
    extent['x_l'] = n[x_var] +  k[x_var] * c_extent[0]
    extent['x_u'] = n[x_var] +  k[x_var] * c_extent[1]
    extent['y_l'] = n[y_var] +  k[y_var] * c_extent[2]
    extent['y_u'] = n[y_var] +  k[y_var] * c_extent[3]
    
    
    #plot the calculated data
    im = ax_im.contourf( Z, cmap = mpl.cm.rainbow, extent = extent.values)
    cbar = fig.colorbar(im)
    cbar.set_label(target_name)
    ax_im.set_xlabel(x_var)
    ax_im.set_ylabel(y_var)
    
    
    
    #for the features that are not ploted we create and set up the sliders
    ax = {}
    sldr = {}
    axcolor = 'lightgoldenrodyellow'
    y_pos = 0.05
    for feature in main_effects:
        ax[feature] = plt.axes([0.15, y_pos, 0.65, 0.03], facecolor=axcolor)
        y_pos = y_pos + 0.05
        valinit = feat_settings[feature]
        sldr[feature] = Slider(ax[feature], feature, -1, 1, valinit = valinit)
 


    def update(val):
        '''
        Update function that is called when slider value is changed.
        '''
        
        #update data
        for feature in main_effects:
            feat_settings[feature] = sldr[feature].val
        #reacalculate feature matrix
        exog = recalculate_feature_matrix(results, feat_settings, x_var, y_var, X, Y, columns)
        
        #predict the values and write data to Z variable
        predicted = results.predict(exog) # evaluation of the function on the grid
        Z = predicted.values.reshape(X.shape)
        
        #prepare for new plot
        ax_im.clear()
        #replot the contour
        im = ax_im.contourf( Z, cmap = mpl.cm.rainbow, extent = extent.values)#, extent = extent)
        ax_im.set_xlabel(x_var)
        ax_im.set_ylabel(y_var)
        
        #update colorbar
        cbar.set_clim(vmax = Z.max(), vmin = Z.min())
        cbar_ticks = np.linspace(round(Z.min()), round(Z.max()), num=10, endpoint=True)
        if Z.max() < 10:
            cbar_ticks = np.linspace(Z.min(), Z.max(), num=10, endpoint=True)
        cbar.set_ticks(cbar_ticks)
        #cbar.update_normal()
        cbar.draw_all()
        fig.canvas.draw()


    #connect the sliders with the change effect
    for feature in main_effects:
        sldr[feature].on_changed(update)
       
    fig.subplots_adjust(left = 0.15, bottom = y_pos + 0.05, top = 0.95, right = 1)
    plt.show()   

def matrix_contourf_predictions():
    pass     

if __name__ == '__main__':

    import chrom_tools
    #select only optimization data
    filter_dict = {'class' : ['optimization'], 'step' : ['MMC']}
    
    #location of files and files
    fname_samples = 'C:/Users/jeromlu2/LukaFiles/04_Programiranje/01_Python/03_IPythonNotebooks/\
temp/IPI/03_DSP Development/04_DSPD/01_Akta_runs/GNT904_data_akta_runs.xlsx'
    fname_orig = 'C:/Users/jeromlu2/LukaFiles/04_Programiranje/01_Python/03_IPythonNotebooks/\
temp/IPI/03_DSP Development/04_DSPD/03_Originator/QTPP data table original.xlsx'
    
    #import of the data
    data = chrom_tools.AnalyticalData(fname_samples, fname_orig, filter_dict)
    originator = data.originator.originator
    
    all_feat_names = ['targ_load','pH','load_cond','buffer_cond','peak_cut_end']
    
    feat_limits = pd.DataFrame({
        'targ_load'   : [50 , 250],
        'pH'         : [6.0, 7.0],
        'load_cond'    : [3, 15],
        'peak_cut_end' : [600, 2000],
        'buffer_cond' : [2, 15]    
    }, index = ['lower', 'upper'])
    
    X_df = data.eluates[all_feat_names]
    
    # transform into matrix of ones and zeros, set max and min value for the transformation
    X_upper_lim = feat_limits.loc['upper']
    X_lower_lim = feat_limits.loc['lower']
    X_df_coded = code_df(X_df, X_upper_lim, X_lower_lim)
    
    #create ful resonse surface
    X_df_coded_RS = create_RS_features(X_df_coded)
    
    target_name = 'eluate_volume'
    
    #Create targets vector
    y_df = data.eluates[target_name]
    
    #scatter_plot_QA(y_df, X_df_coded)
    results, model, history, cross_scores, Q2_scores = stepwise_feat_selection(y_df, X_df_coded_RS, restrict= True, cv = 10)
    #plot_coeff_effects(results, 'pvalues')
    plot_history_of_stepwise(target_name, history, cross_scores, Q2_scores)
    
    columns = X_df.columns
    #prediction_profiler(results, feat_limits, columns)
    
    x_var = 'pH'
    y_var = 'load_cond'
    c_extent = [-1, 1, -1, 1]
    
    prediction_contour_plot(results, x_var, y_var, feat_limits, c_extent, columns)
    
    #fig, ax_list, sldrs = prediction_profiler(results, feat_limits, columns);
    
 