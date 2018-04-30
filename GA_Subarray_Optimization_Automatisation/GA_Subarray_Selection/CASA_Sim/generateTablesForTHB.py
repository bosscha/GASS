import os
import glob
import numpy

configs = ['ACA-std.cfg']
for i in sorted(glob.glob('C43*cfg')):
    configs.append(i)

weightings = ['natural', 'briggs', 'uniform']

# for config in configs:
#     for weighting in weightings:
#     
#         os.system('python arrayConfiguration.py -i '+config+' -t casa -l 1h -a -0.5h -d -23 -w '+weighting+' -j report_'+config.replace(".cfg", "")+'_'+weighting+'.json')

reports = glob.glob('report_*')

data = {}
configs1 = []

for i in range(len(reports)):

    f = open(reports[i])
    fc = eval(f.read())
    f.close()

    ij = fc['inputCfg'].replace('.cfg', '')+'_'+fc['weighting']
    data[ij] = fc

    ij = tuple([fc['inputCfg'], fc['maximumBaseline']])
    if ij not in configs1:
        configs1.append(ij)

configs1 = [j[0] for j in sorted(configs1, key=lambda x:x[1])]

configs2 = [configs1[0:6], configs1[7:]]

for configs in configs2:

    line = '\\begin{tabular}{|'
    for i in range(len(configs)+1):
        line += 'l|'
    line += '}\n\\hline\nConfiguration'
    for i in range(len(configs)):
        config1 = configs[i].replace('.cfg', '')
        if config1 == 'ACA-std':
            config1 = '7 m'
        line += ' & ' + config1
    line += ' \\\\\n\hline\nMinimum baseline (m)'
    for i in range(len(configs)):
        ij = configs[i].replace('.cfg', '')+'_briggs'
        line += ' & ' + str(round(data[ij]['minimumBaseline'], 1))
    line += ' \\\\\n5th percentile or L05 (m)'
    for i in range(len(configs)):
        ij = configs[i].replace('.cfg', '')+'_briggs'
        line += ' & ' + str(round(data[ij]['percentileBaseline']['L05'], 1))
    line += ' \\\\\n80th percentile or L80 (m)'
    for i in range(len(configs)):
        ij = configs[i].replace('.cfg', '')+'_briggs'
        line += ' & ' + str(round(data[ij]['percentileBaseline']['L80'], 1))
    line += ' \\\\\nMaximum baseline (m)'
    for i in range(len(configs)):
        ij = configs[i].replace('.cfg', '')+'_briggs'
        line += ' & ' + str(round(data[ij]['maximumBaseline'], 1))
    line += ' \\\\\n\hline\n\\end{tabular}\n'

    print line

print ''

for configs in configs2:

    line = '\\begin{tabular}{|'
    for i in range(len(configs)+1):
        line += 'l|'
    line += '}\n\\hline\nConfiguration'
    for i in range(len(configs)):
        config1 = configs[i].replace('.cfg', '')
        if config1 == 'ACA-std':
            config1 = '7 m'
        line += ' & ' + config1

    for weighting in weightings:
        weighting1 = weighting
        if weighting == 'briggs':
            weighting1 += ' ($R$ = 0.5)'
        line += ' \\\\\n\hline\n'+weighting1+' '
        for i in range(len(configs)):
            ij = configs[i].replace('.cfg', '')+'_'+weighting
            line += ' & ' + str(round(data[ij]['sidelobeLevel'], 1))+'\%'

    line += ' \\\\\n\hline\n\\end{tabular}\n'

    print line
