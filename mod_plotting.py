__author__ = 'boris'
from scipy.stats.mstats import winsorize
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy
import math
import mod_lib
import mod_utils
import operator
import os
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines


def plot_mutated_nts_pie(libraries, out_prefix, subtract_background=False, subtract_control=False, exclude_constitutive=False):
    #Makes an array of pie charts, 1 per library
    if subtract_background or subtract_control:
        #if subtracting background, need to only look at those which have a defined control
        libraries = [library for library in libraries if library.lib_settings.sample_name in
                     library.experiment_settings.get_property('experimentals')]
    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mutated_nts_count = library.count_mutation_rates_by_nucleotide(subtract_background=subtract_background, subtract_control=subtract_control,
                                                                       exclude_constitutive=exclude_constitutive)
        labels = sorted(mutated_nts_count.keys())
        sizes = numpy.array([mutated_nts_count[nt] for nt in labels])
        total = float(sum(sizes))
        sizes = sizes/total
        merged_labels = ['%s %.3f' % (labels[i], sizes[i]) for i in range(len(sizes))]
        plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
        plot.set_title(library.lib_settings.sample_name)
        plot_index += 1
    if subtract_background:
        plt.suptitle('background-subtracted mutation rate fractions')
    if subtract_control:
        plt.suptitle('control-subtracted mutation rate fractions')
    else:
        plt.suptitle('mutation rate fractions')
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_mutation_breakdown_pie(libraries, out_prefix, exclude_constitutive=False):
    #Makes an array of pie charts, 4 pers library, with types of mutations for each nt
    num_subplots = len(libraries)*4
    num_plots_wide = 4
    num_plots_high = len(libraries)
    fig = plt.figure(figsize=(16, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        mutated_nts_count = library.count_mutation_types_by_nucleotide(exclude_constitutive=exclude_constitutive)
        for nt in 'ATCG':
            plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
            sorted_muts = sorted(mutated_nts_count[nt].items(), key=operator.itemgetter(1), reverse=True)
            labels = [pair[0] for pair in sorted_muts[:4]]
            others = [pair[0] for pair in sorted_muts[4:]]
            sizes = [mutated_nts_count[nt][mut_type] for mut_type in labels]
            labels.append('other')
            other_sum = sum([mutated_nts_count[nt][mut_type] for mut_type in others])
            sizes.append(other_sum)
            sizes = numpy.array(sizes)
            total = float(sum(sizes))
            sizes = sizes/total
            merged_labels = ['%s %.3f' % (labels[i], sizes[i]) for i in range(len(sizes))]
            plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
            plot.set_title(library.lib_settings.sample_name)
            plot_index += 1
    plt.suptitle('mutation rate type fractions')
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_mutation_rate_cdfs(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False):
    #Makes 2 CDF plots. One of all libraries, showing the coverage-normalized mutation rates
    # and one showing background-subtracted mutation rates

    fig = plt.figure(figsize=(24, 16))
    plots = []
    plot = fig.add_subplot(231)
    plots.append(plot)
    colormap = plt.get_cmap('spectral')
    colorindex = 0
    for library in libraries:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("mutation rate")
    plot.set_title('raw mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(232)
    plots.append(plot)
    colorindex = 0  #Set this to 1 instead to keep the color scheme the same after No DMS has been excluded.
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=True, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("background-subtracted mutation rate")
    plot.set_title('normalized mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(233)
    plots.append(plot)
    colorindex = 0  #Set this to 1 instead to keep the color scheme the same after No DMS has been excluded. 
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=True,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("control-subtracted mutation rate")
    plot.set_title('control normalized mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(234)
    plots.append(plot)
    colormap = plt.get_cmap('spectral')
    colorindex = 0
    for library in libraries:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("log10 mutation rate")
    plot.set_title('raw mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    plot = fig.add_subplot(235)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=True, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("background-subtracted log10 mutation rate")
    plot.set_title('normalized mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    plot = fig.add_subplot(236)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=True,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("control-subtracted log10 mutation rate")
    plot.set_title('control normalized mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    for plot in plots:
        plot.set_ylabel("cumulative fraction of %s nucleotides" % (nucleotides_to_count))
        plot.set_ylim(0, 1)
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_changes_vs_control_interactive(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
                                        max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of magnitude of change (experimental-control)
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    from bokeh.plotting import figure, output_file, show, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict

    # output to static HTML file
    output_file("%s.html" % (out_prefix))
    plot_figs=[]

    for library in libraries:
        mag_change, fold_change, annotation = [], [], []
        prot_mag_change, prot_fold_change, prot_annotation = [], [], []
        deprot_mag_change, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        #print 'deprotected ', nucleotide
                        deprot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        #print 'protected ', nucleotide
                        prot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        source = ColumnDataSource(data=dict(x = mag_change, y = fold_change, label = annotation))
        prot_source = ColumnDataSource(data=dict(x = prot_mag_change, y = prot_fold_change, label = prot_annotation))
        deprot_source = ColumnDataSource(data=dict(x = deprot_mag_change, y = deprot_fold_change,
                                                   label = deprot_annotation))
        TOOLS = "pan,wheel_zoom,reset,save,hover"
        PlotFig = figure(x_axis_label = "mutation rate [%s] - [%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_label = "fold change [%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_type="log", tools=TOOLS, toolbar_location="right")
        PlotFig.circle("x", "y", size = 5, source=source, color=mod_utils.bokeh_black)
        PlotFig.circle("x", "y", size = 5, source=prot_source, color=mod_utils.bokeh_vermillion)
        PlotFig.circle("x", "y", size = 5, source=deprot_source, color=mod_utils.bokeh_bluishGreen)
        PlotFig.x_range = Range1d(start=-0.2, end=0.2)
        PlotFig.y_range = Range1d(start=.001, end=100)

        #adjust what information you get when you hover over it
        Hover = PlotFig.select(dict(type=HoverTool))
        Hover.tooltips = OrderedDict([("nuc", "@label")])
        plot_figs.append([PlotFig])
    p = gridplot(plot_figs)
    show(p)

def ma_plots_interactive(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
                         max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of average magnitude of signal (experimental+control)/2
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    from bokeh.plotting import figure, output_file, show, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict

    # output to static HTML file
    output_file("%s.html" % (out_prefix))
    plot_figs=[]

    for library in libraries:
        mag, fold_change, annotation = [], [], []
        prot_mag, prot_fold_change, prot_annotation = [], [], []
        deprot_mag, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    avg_mutation_rate = (nucleotide.mutation_rate+nucleotide.get_control_nucleotide().mutation_rate)/2.0
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag.append(avg_mutation_rate)
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag.append(avg_mutation_rate)
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag.append(avg_mutation_rate)
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        source = ColumnDataSource(data=dict(x = mag, y = fold_change, label = annotation))
        prot_source = ColumnDataSource(data=dict(x = prot_mag, y = prot_fold_change, label = prot_annotation))
        deprot_source = ColumnDataSource(data=dict(x = deprot_mag, y = deprot_fold_change,
                                                   label = deprot_annotation))
        TOOLS = "pan,wheel_zoom,reset,save,hover"
        PlotFig = figure(x_axis_label = "avg signal ([%s] + [%s])/2" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_label = "fold change [%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_type="log", x_axis_type="log", tools=TOOLS, toolbar_location="right")
        PlotFig.circle("x", "y", size = 5, source=source, color=mod_utils.bokeh_black)
        PlotFig.circle("x", "y", size = 5, source=prot_source, color=mod_utils.bokeh_vermillion)
        PlotFig.circle("x", "y", size = 5, source=deprot_source, color=mod_utils.bokeh_bluishGreen)
        PlotFig.x_range = Range1d(start=0.00001, end=1)
        PlotFig.y_range = Range1d(start=.001, end=100)

        #adjust what information you get when you hover over it
        Hover = PlotFig.select(dict(type=HoverTool))
        Hover.tooltips = OrderedDict([("nuc", "@label")])
        plot_figs.append([PlotFig])
    p = gridplot(plot_figs)
    show(p)

def plot_changes_vs_control(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
                            max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library make a plot of magnitude of change (experimental-control)
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
    """
    output_file = "%s.pdf" % (out_prefix)
    plot_figs=[]

    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mag_change, fold_change, annotation = [], [], []
        prot_mag_change, prot_fold_change, prot_annotation = [], [], []
        deprot_mag_change, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        plot.set_xlabel("[%s] - [%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_ylabel("[%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),  fontsize = 8)
        plot.set_yscale('log')
        plot.scatter(mag_change, fold_change, color=mod_utils.black, s=2)
        plot.scatter(prot_mag_change, prot_fold_change, color=mod_utils.vermillion, s=2)
        plot.scatter(deprot_mag_change, deprot_fold_change, color=mod_utils.bluishGreen, s=2)
        plot.set_xlim(-0.2,0.2)
        plot.set_ylim(.001,100)

        plot_figs.append(plot)
        plot_index+=1
    plt.savefig(output_file, transparent='True', format='pdf')

def ma_plots(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
             max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of magnitude of signal (experimental+control)/2
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    output_file = "%s.pdf" % (out_prefix)
    plot_figs=[]

    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mag, fold_change, annotation = [], [], []
        prot_mag, prot_fold_change, prot_annotation = [], [], []
        deprot_mag, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    avg_mutation_rate = (nucleotide.mutation_rate+nucleotide.get_control_nucleotide().mutation_rate)/2.0
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag.append(avg_mutation_rate)
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag.append(avg_mutation_rate)
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag.append(avg_mutation_rate)
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        plot.set_xlabel("([%s] + [%s])/2" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_ylabel("[%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_yscale('log')
        plot.set_xscale('log')
        plot.scatter(mag, fold_change, color=mod_utils.black, s=3)
        plot.scatter(prot_mag, prot_fold_change, color=mod_utils.vermillion, s=5)
        plot.scatter(deprot_mag, deprot_fold_change, color=mod_utils.bluishGreen, s=5)
        plot.set_xlim(0.00001,1)
        plot.set_ylim(.001,100)
        plot_figs.append(plot)
        plot_index+=1
    plt.savefig(output_file, transparent='True', format='pdf')

def highlight_structure(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of magnitude of signal (experimental+control)/2
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """

    for library in libraries:

        protected_nucleotides = library.get_changed_nucleotides('protected', confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
        num_protected = 0
        for rRNA in protected_nucleotides:
                num_protected += len(protected_nucleotides[rRNA])
        deprotected_nucleotides = library.get_changed_nucleotides('deprotected', confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
        num_deprotected = 0
        for rRNA in deprotected_nucleotides:
            num_deprotected += len(deprotected_nucleotides[rRNA])

        if num_protected>0 or num_deprotected>0:
            output_file = open(os.path.join(out_prefix, "%s.txt" % (library.lib_settings.sample_name)), 'w')
            reference_pymol_script_file = open(library.experiment_settings.get_property('pymol_base_script'), 'rU')
            for line in reference_pymol_script_file:
                if line.startswith('#<insert nucleotide highlighting here>'):
                    if num_protected>0:
                        rRNA_selections = []
                        for rRNA in protected_nucleotides:
                            if len(protected_nucleotides[rRNA])>0:
                                rRNA_selections.append('%s and resi %s' % (rRNA, '+'.join([str(nucleotide.position) for
                                                                                         nucleotide in protected_nucleotides[rRNA]])))
                        outline = 'create protected_nucleotides, %s\n' % (' or '.join(rRNA_selections))
                        output_file.write(outline)
                    if num_deprotected>0:
                        rRNA_selections = []
                        for rRNA in deprotected_nucleotides:
                            if len(deprotected_nucleotides[rRNA])>0:
                                rRNA_selections.append('%s and resi %s' % (rRNA, '+'.join([str(nucleotide.position) for
                                                                                         nucleotide in deprotected_nucleotides[rRNA]])))
                        outline = 'create deprotected_nucleotides, %s\n' % (' or '.join(rRNA_selections))
                        output_file.write(outline)
                elif line.startswith('#<color groups here>'):
                    if num_protected>0:
                        output_file.write('color vermillion, protected_nucleotides\n')
                    if num_deprotected>0:
                        output_file.write('color bluish_green, deprotected_nucleotides\n')
                elif line.startswith('#<show spheres for changing nucleotides here>'):
                    if num_protected>0:
                        output_file.write('show spheres, protected_nucleotides\n')
                    if num_deprotected>0:
                        output_file.write('deprotected_nucleotides\n')
                else:
                    output_file.write(line)
            reference_pymol_script_file.close()
            output_file.close()


def color_by_change(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False, subtract_background=False):
    for library in libraries:
        log_fold_changes = {}
        maxval = 0.0
        minval = 0.0
        for rRNA_name in library.rRNA_mutation_data:
            log_fold_changes[rRNA_name] = {}
            for nucleotide in library.rRNA_mutation_data[rRNA_name].nucleotides:
                if library.rRNA_mutation_data[rRNA_name].nucleotides[nucleotide].identity in nucleotides_to_count and \
                        library.rRNA_mutation_data[rRNA_name].nucleotides[nucleotide].\
                                get_control_fold_change_in_mutation_rate(subtract_background=subtract_background) not in [0.0, float('inf'), float('-inf')]:
                    if exclude_constitutive and library.rRNA_mutation_data[rRNA_name].nucleotides[nucleotide].exclude_constitutive:
                        pass
                    else:
                        log_fold_changes[rRNA_name][nucleotide] = math.log(library.rRNA_mutation_data[rRNA_name].nucleotides[nucleotide].get_control_fold_change_in_mutation_rate(subtract_background=subtract_background), 10)
                        maxval = max(maxval, log_fold_changes[rRNA_name][nucleotide])
                        minval = min(minval, log_fold_changes[rRNA_name][nucleotide])
                else:
                    pass

        absmax = max(abs(maxval), abs(minval))
        output_file = open(os.path.join(out_prefix, "%s.txt" % (library.lib_settings.sample_name)), 'w')
        reference_pymol_script_file = open(library.experiment_settings.get_property('pymol_base_script_colorchange'), 'rU')
        for line in reference_pymol_script_file:
            if line.startswith('#<insert b-factors>'):
                output_file.write('python\n')
                output_file.write('cmd.alter(\'all\', \'b=0.0\')\n')
                for rRNA_name in log_fold_changes:
                    for nucleotide in log_fold_changes[rRNA_name]:
                        output_file.write('cmd.alter(\''+rRNA_name+' and resi '+str(nucleotide)+'\', \'b=float("'+str(log_fold_changes[rRNA_name][nucleotide])+'")\')\n')
                output_file.write('python end\n')
            elif line.startswith('#<insert spectrum>'):
                output_file.write('spectrum b, bluish_green white vermillion, minimum='+str(-absmax)+', maximum='+str(absmax)+'\n')
                output_file.write('ramp_new scale, E.c.16S_rRNA, ['+str(-absmax)+',0,'+str(absmax)+'], [bluish_green, white, vermillion]')
                output_file.write('ramp_new scale, E.c.23S__rRNA, ['+str(-absmax)+',0,'+str(absmax)+'], [bluish_green, white, vermillion]')
            else:
                output_file.write(line)
        reference_pymol_script_file.close()
        output_file.close()


        output_file = open(os.path.join(out_prefix, "%s_xtreme.txt" % (library.lib_settings.sample_name)), 'w')
        reference_pymol_script_file = open(library.experiment_settings.get_property('pymol_base_script_colorchange'), 'rU')
        for line in reference_pymol_script_file:
            if line.startswith('#<insert b-factors>'):
                output_file.write('python\n')
                output_file.write('cmd.alter(\'all\', \'b=0.0\')\n')
                for rRNA_name in log_fold_changes:
                    for nucleotide in log_fold_changes[rRNA_name]:
                        if log_fold_changes[rRNA_name][nucleotide] > 0.699 or log_fold_changes[rRNA_name][nucleotide] < -0.177: #previously 0.477 and -0.398
                            output_file.write('cmd.alter(\''+rRNA_name+' and resi '+str(nucleotide)+'\', \'b=float("'+str(log_fold_changes[rRNA_name][nucleotide])+'")\')\n')
                        else:
                            pass
                output_file.write('python end\n')
            elif line.startswith('#<insert spectrum>'):
                output_file.write('spectrum b, bluish_green white vermillion, minimum='+str(-absmax)+', maximum='+str(absmax)+'\n')
                output_file.write('ramp_new scale, E.c.16S_rRNA, ['+str(-absmax)+',0,'+str(absmax)+'], [bluish_green, white, vermillion]')
                output_file.write('ramp_new scale, E.c.23S__rRNA, ['+str(-absmax)+',0,'+str(absmax)+'], [bluish_green, white, vermillion]')
            else:
                output_file.write(line)
        reference_pymol_script_file.close()
        output_file.close()


        for rRNA_name in log_fold_changes:
            xrna_file = open(os.path.join(out_prefix, "%s_%s.xrna" % (library.lib_settings.sample_name, rRNA_name[4:])), 'w')
#            reference_xrna_script_file = open(library.experiment_settings.get_property('xrna_base_script_%s' % (rRNA_name[4:])), 'rU')
            reference_xrna_script_file = open(library.experiment_settings.get_property('xrna_base_script_colorchange'), 'rU')
            for line in reference_xrna_script_file:
                if line.startswith('#<insert colored nucleotides>'):
                    for nucleotide in log_fold_changes[rRNA_name]:
                        if log_fold_changes[rRNA_name][nucleotide] > 0.699:
                            xrna_file.write('<Nuc RefID=\''+str(nucleotide)+'\' Color=\'ff'+format(int(255*(1-((log_fold_changes[rRNA_name][nucleotide])/absmax))), '02x')*2+'\'>\n')
                            xrna_file.write('<NucSymbol>\n')
                            xrna_file.write('a 0.0 0.0 5.0 0.0 360.0 0.5 0  0\n')
                            xrna_file.write('</NucSymbol>\n')
                            xrna_file.write('</Nuc>\n')
                        elif log_fold_changes[rRNA_name][nucleotide] < -0.177:
                            xrna_file.write('<Nuc RefID=\''+str(nucleotide)+'\' Color=\''+format(int(255*(1-abs((log_fold_changes[rRNA_name][nucleotide])/absmax))), '02x')*2+'ff\'>\n')
                            xrna_file.write('<NucSymbol>\n')
                            xrna_file.write('a 0.0 0.0 5.0 0.0 360.0 0.5 0  0\n')
                            xrna_file.write('</NucSymbol>\n')
                            xrna_file.write('</Nuc>\n')
                        else:
                            pass
#                            xrna_file.write('<Nuc RefID=\''+str(nucleotide)+'\' Color=\'ffffff\'>\n')
#                            xrna_file.write('<NucSymbol>\n')
#                            xrna_file.write('a 0.0 0.0 5.0 0.0 360.0 0.5 0  0\n')
#                            xrna_file.write('</NucSymbol>\n')
#                            xrna_file.write('</Nuc>\n')
                else:
                    xrna_file.write(line)
            reference_xrna_script_file.close()
            xrna_file.close()


def generate_roc_curves(tp_tn_annotations, genome_fasta, outprefix, libraries, rRNA, nucs_to_count):
    def winsorize_norm_chromosome_data(mut_density, chromosome, genome_dict, nucs_to_count, to_winsorize = False, low = 0, high = 0.95):
        """


        :param read_5p_ends:
        :param chromosome:
        :param strand:
        :param genome_dict:
        :param nucs_to_count:
        :param low:
        :param high:
        :return: an array (now zero-indexed from 1-indexed) of densities for the given chromosome on the given strand, winsorized, and only for the given nucleotides
        """
        max_position = max(mut_density[chromosome].nucleotides.keys())
        density_array =numpy.array([0.0] * max_position)
        for position in mut_density[chromosome].nucleotides.keys():
            if genome_dict[chromosome][position-1] in nucs_to_count:
                density_array[position-1] = mut_density[chromosome].nucleotides[position].mutation_rate
        if to_winsorize:
            winsorize(density_array, limits = (low, 1-high), inplace = True)
        normed_array = density_array/float(max(density_array))
        return  normed_array

    def get_tp_tn(tp_tn_file):
        TP = set()
        TN = set()
        f = open(tp_tn_file)
        for line in f:
            ll= line.strip('\n').split('\t')
            if ll[2] == 'TP':
                TP.add(int(ll[0]))
            if ll[2] =='TN':
                TN.add(int(ll[0]))
        f.close()
        return TP, TN

    def call_positives(density_array, chromosome, genome_dict, nucs_to_count, cutoff):
        """

        :param density_array:
        :return:a set of called positive positions
                I've reverted these to 1-indexed to match the TP and TN calls from the structures
        """
        positives = set()

        for i in range(len(density_array)):
            if genome_dict[chromosome][i] in nucs_to_count:
                if density_array[i] >= cutoff:
                    positives.add(i+1)#adding 1 not necessary for RT stops, since the modified nucleotide is the one 1 upstream of the RT stop!!!

        return positives

    def plot_ROC_curves(roc_curves, title, out_prefix):
        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)#first a pie chart of mutated nts
        colormap = plt.get_cmap('spectral')
        color_index = 0
        for name in sorted(roc_curves.keys()):
            x, y = roc_curves[name]
            area_under_curve = numpy.trapz(numpy.array(y[::-1])/100., x=numpy.array(x[::-1])/100.)
            plot.plot(x, y, lw =2, label = '%s   %.3f' % (name, area_under_curve), color = colormap(color_index/float(len(roc_curves))))
            color_index +=1
        plot.plot(numpy.arange(0,100,0.1), numpy.arange(0,100,0.1), lw =1, ls = 'dashed', color = mod_utils.black, label = 'y=x')
        plot.set_xlabel('False positive rate (%) (100-specificity)')
        plot.set_ylabel('True positive rate (%) (sensitivity)')
        plot.set_title(title)
        lg=plt.legend(loc=4,prop={'size':10}, labelspacing=0.2)
        lg.draw_frame(False)
        plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
        plt.clf()

    sample_names = [library.lib_settings.sample_name for library in libraries]
    mutation_densities = [library.rRNA_mutation_data for library in libraries]


    genome_dict = genome_fasta
    normed_density_arrays = [winsorize_norm_chromosome_data(mutation_density, rRNA, genome_dict, nucs_to_count) for mutation_density in mutation_densities]
    real_tp, real_tn = get_tp_tn(tp_tn_annotations)
    roc_curves = {}
    for sample_name in sample_names:
        roc_curves[sample_name] = [[],[]]#x and y value arrays for each

    stepsize = 0.0001
    for cutoff in numpy.arange(0,1.+5*stepsize, stepsize):
        for i in range(len(sample_names)):
            called_p = call_positives(normed_density_arrays[i], rRNA, genome_dict, nucs_to_count, cutoff)
            num_tp_called = len(called_p.intersection(real_tp))#how many true positives called at this cutoff
            num_fp_called = len(called_p.intersection(real_tn))#how many fp positives called at this cutoff
            roc_curves[sample_names[i]][1].append(100.*num_tp_called/float(len(real_tp)))#TP rate on y axis
            roc_curves[sample_names[i]][0].append(100.*num_fp_called/float(len(real_tn)))#FP rate on x axis

    plot_ROC_curves(roc_curves, rRNA, outprefix)

def parse_functional_groups(groups_file, delimiter='\t'):
    return_dict = defaultdict(list)
    f = open(groups_file, 'rU')
    lines  = f.readlines()
    headers = lines[0].strip('\n').split(delimiter)
    for line in lines[1:]:
        ll= line.strip('\n').split(delimiter)
        for i in range(len(ll)):
            if not ll[i].strip() =='':
                return_dict[headers[i]].append(ll[i])
    f.close()
    return return_dict

def plot_functional_group_changes(libraries, out_prefix, groups_file, nucleotides_to_count='ATCG', exclude_constitutive=False,
                            max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library make a plot of log10 fold change (experimental/control) between different functional groups.
            both as a CDF and as a violin plot
    """
    functional_groups = parse_functional_groups(groups_file)
    group_names = sorted(functional_groups.keys())
    for library in libraries:
        cdf_file = "%s_%s_CDF.pdf" % (out_prefix, library.lib_settings.sample_name)
        num_plots_wide = 1
        num_plots_high = 1
        fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
        plot_index =1
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        colorindex = 0
        all_fold_changes = library.list_fold_changes(nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)
        hist, bin_edges = numpy.histogram(all_fold_changes, bins=10000)
        cum_hist = numpy.cumsum(hist)
        cum_hist = cum_hist/float(max(cum_hist))
        plot.plot(bin_edges[:-1], cum_hist, color=mod_utils.colors[0], label='all %d (K-S P)' % (len(all_fold_changes)), lw=2)
        for group_name in group_names:
            colorindex+=1
            group_fold_changes = [nucleotide.get_control_fold_change_in_mutation_rate() for nucleotide in
                                  library.get_nucleotides_from_list(functional_groups[group_name],
                                                                   nucleotides_to_count=nucleotides_to_count,
                                                                   exclude_constitutive=exclude_constitutive)]
            d, p = stats.ks_2samp(all_fold_changes, group_fold_changes)
            hist, bin_edges = numpy.histogram(group_fold_changes, bins=10000)
            cum_hist = numpy.cumsum(hist)
            cum_hist = cum_hist/float(max(cum_hist))
            #print group_fold_changes
            plot.plot(bin_edges[:-1], cum_hist, color=mod_utils.colors[colorindex],label='%s %d (%f)' % (group_name, len(group_fold_changes), p), lw=2)

        lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
        lg.draw_frame(False)
        plot.set_ylabel("cumulative nucleotide fraction", fontsize = 8)
        plot.set_xlabel("[%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),  fontsize = 8)
        plot.set_xscale('log')
        plot.set_xlim(.1,10)
        plot.set_ylim(0, 1)

        plt.savefig(cdf_file, transparent='True', format='pdf')
