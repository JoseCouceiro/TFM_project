import os
import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import streamlit as st

class Graphs():
    """
    This class contains functions necessary to plot a classification report and a confussion matrix.
    The classification report functions were written by HYRY as credited, and only slighly modified.
    They need the output of 'metrics.sklearn.classification_report' as a text file to work.
    """

    def __init__(self):
        with open(os.path.join('..', '..', 'data', '03_primary', 'code_to_label_dic.json')) as codes:
            self.label_dic = json.load(codes)        

    def show_values(self, pc, fmt="%.2f", **kw):
        '''
        Heatmap with text in each cell with matplotlib's pyplot
        Source: https://stackoverflow.com/a/25074150/395857 
        By HYRY
        '''
        pc.update_scalarmappable()
        ax = pc.axes
        for p, color, value in zip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
            x, y = p.vertices[:-2, :].mean(0)
            if np.all(color[:3] > 0.5):
                color = (0.0, 0.0, 0.0)
            else:
                color = (1.0, 1.0, 1.0)
            ax.text(x, y, fmt % value, ha="center", va="center", color=color, **kw)

    def cm2inch(self, *tupl):
        '''
        Specify figure size in centimeter in matplotlib
        Source: https://stackoverflow.com/a/22787457/395857
        By gns-ank
        '''
        inch = 2.54
        if type(tupl[0]) == tuple:
            return tuple(i/inch for i in tupl[0])
        else:
            return tuple(i/inch for i in tupl)
        
    def heatmap(self, AUC, title, xlabel, ylabel, xticklabels, yticklabels, figure_width=40, figure_height=20, correct_orientation=False, cmap='RdBu'):
        '''
        Inspired by:
        - https://stackoverflow.com/a/16124677/395857 
        - https://stackoverflow.com/a/25074150/395857
        '''

        # Plot it out
        fig, ax = plt.subplots()    
        #c = ax.pcolor(AUC, edgecolors='k', linestyle= 'dashed', linewidths=0.2, cmap='RdBu', vmin=0.0, vmax=1.0)
        c = ax.pcolor(AUC, edgecolors='k', linestyle= 'dashed', linewidths=0.2, cmap=cmap)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(AUC.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(AUC.shape[1]) + 0.5, minor=False)

        # set tick labels
        #ax.set_xticklabels(np.arange(1,AUC.shape[1]+1), minor=False)
        ax.set_xticklabels(xticklabels, minor=False)
        ax.set_yticklabels(yticklabels, minor=False)

        # set title and x/y labels
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)      

        # Remove last blank column
        plt.xlim( (0, AUC.shape[1]) )

        # Turn off all the ticks
        ax = plt.gca()    
        for t in ax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False

        # Add color bar
        plt.colorbar(c)
        # Add text in each cell 
        self.show_values(c)

        # Proper orientation (origin at the top left instead of bottom left)
        if correct_orientation:
            ax.invert_yaxis()
            ax.xaxis.tick_top()       

        # resize 
        fig = plt.gcf()
        fig.set_size_inches(self.cm2inch(figure_width, figure_height))

    def plot_classification_report(self, classification_report, title='Classification report ', cmap='RdBu'):
        '''
        Plot scikit-learn classification report.
        Extension based on https://stackoverflow.com/a/31689645/395857 
        '''
        lines = classification_report.split('\n')

        #classes = []
        plotMat = []
        support = []
        class_names = []
        for line in lines[2 : (len(lines))-5]:
            t = line.strip().split()
            if len(t) < 2: continue
            v = [float(x) for x in t[1: len(t) - 1]]
            support.append(int(t[-1]))
            class_names.append(t[0])
            plotMat.append(v)
        
        class_names = [self.label_dic[str(key)] for key in class_names]
        
        xlabel = 'Metrics'
        ylabel = 'Classes'
        xticklabels = ['Precision', 'Recall', 'F1-score']
        yticklabels = ['{0} ({1})'.format(class_names[idx], sup) for idx, sup  in enumerate(support)]
        figure_width = 30
        figure_height = len(class_names) + 8.4
        correct_orientation = False
        self.heatmap(np.array(plotMat), title, xlabel, ylabel, xticklabels, yticklabels, figure_width, figure_height, correct_orientation, cmap=cmap)

    def plot_confusion_matrix(self, y_pred, y_test):
        """ Function that calculates a confussion matrix and plots it into a graphic.
        Args: predictions, true_values.
        Output: a matplotlib figure.
        """

        try:
            #Gather the total number of labels
            labels = list(set(y_pred).union(set(np.unique(y_test))))
            y_labels = [self.label_dic[str(y)].split()[0] for y in labels]

            #Plot
            fig, ax = plt.subplots(figsize=(10,10))
            disp=metrics.ConfusionMatrixDisplay.from_predictions(
                y_test,y_pred, normalize='true', display_labels=y_labels, xticks_rotation='vertical', values_format='.1f', ax=ax)
            
            return disp.figure_
        except TypeError:
            st.error("Something went wrong! Please, make sure there are no 'nan' in the label column")
