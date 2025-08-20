import numpy as np
import jpype
import os
from sympy import Symbol, simplify
import re

def setFuntions(selectedFunctions):
    nameFileNodeConfig = "NodeConfig.txt"
    text_file = open(nameFileNodeConfig, "r")
    lines = text_file.readlines()
    text_file.close()
    idx = [i for i, ss in enumerate(lines) if "realFunctions" in ss][0]
    lineTemp = 'realFunctions = '
    for nameFuntion in selectedFunctions:
        lineTemp += nameFuntion + ', '
    lines[idx] = lineTemp[:-2] + '\n'
    text_file = open(nameFileNodeConfig, "w")
    lines = text_file.write(''.join(lines))
    text_file.close()


def randrange(nsamples, lims):
    a, b = lims
    return (b - a) * np.random.rand(*nsamples) + a


def py2java(x, trainData=False):  # matrix [nxm]python to [mxn]java
    if len(x.shape) == 1:
        if trainData:
            res = jpype.JArray(jpype.JDouble, 2)([x.tolist()])
        else:
            res = jpype.JArray(jpype.JDouble, 1)(x.tolist())
    elif len(x.shape) == 2:
        if np.min(x.shape) == 1:
            if trainData:
                res = jpype.JArray(jpype.JDouble, 2)([x.flatten().tolist()])
            else:
                res = jpype.JArray(jpype.JDouble, 1)(x.flatten().tolist())
        else:
            res = jpype.JArray(jpype.JDouble, 2)(x.T.tolist())
    return res


def java2py(x):
    return np.asarray(x, dtype=np.float)


def getTree2plot(t, count=-1, ndecimal=3):
    def roundNumbers(p):
        if hasattr(p, "constant"):
            temp = '%.' + '%df' % (ndecimal)
            y = temp % (p.constant)
        else:
            y = p.name()
        return y

    countOld = count
    count = count + 1
    output = ''
    if count == 0:
        t = t.getKid(0)
        fill = "#136ed4"  # color nodo
        output = "digraph program {\nnode [style=filled]"
        if t.nKids()>0:
            output += ('%d [label="%s", fillcolor="%s"] ;\n' % (count, t.name(), fill))
            tempOutput, _ = getTree2plot(t, count)  # root node
        else:
            fill = "#60a6f6"
            output += ('%d [label="%s", fillcolor="%s"] ;\n' % (count, roundNumbers(t), fill))
            tempOutput = ''
        output += tempOutput
        output += "}"
    else:
        if t.nKids() > 0:  # is node
            output = ''
            for i in range(t.nKids()):
                if t.getKid(i).nKids() > 0:
                    fill = "#136ed4"  # color nodo
                    output += '%d [label="%s", fillcolor="%s"] ;\n' % (count, t.getKid(i).name(), fill)
                    tempOutput, counttemp = getTree2plot(t.getKid(i), count)
                    output += tempOutput
                    output += '%d -> %d ;\n' % (countOld, count)
                    count = counttemp
                else:
                    fill = "#60a6f6"
                    output += '%d [label="%s", fillcolor="%s"] ;\n' % (count, roundNumbers(t.getKid(i)), fill)
                    output += '%d -> %d ;\n' % (countOld, count)
                    count = count + 1

    if countOld == -1:
        return output
    else:
        return output, count

def getReducedAlgebraic(t):
    ss = getAlgebraic(t)
    symbolsfound = list(set(re.findall('X\d*',ss)))
    listsymbols = list()
    for s1 in symbolsfound:
        listsymbols.append([s1, 'Symbol(\''+s1+'\')'])
    vars().update(dict(listsymbols))
    return str(simplify(ss))

def getAlgebraic(t, count=-1):
    countOld = count
    count = count + 1
    
    if count == 0:
        output = ''
        t1 = t.getKid(0)
        if t1.nKids()>0:
            tempOutput, _ = getAlgebraic(t1, count)  # root node
            if len(tempOutput) == 1:
                output = t1.name() + '(' + tempOutput[0] + ')'
            elif len(tempOutput) == 2:
                output = '( ' + tempOutput[0] + ' ' + t1.name() + ' ' + tempOutput[1] + ' )'
            else:
                print('Error. Only binary treeccs')
        else:
            output = t1.name()
    else:
        output = list()
        for i in range(t.nKids()):
            if t.getKid(i).nKids() > 0:
                tempOutput, counttemp = getAlgebraic(t.getKid(i), count)
                if len(tempOutput) == 1:
                    if t.getKid(i).name() == 'sq':
                        temp = tempOutput[0] + '**2'
                    else:
                        temp = t.getKid(i).name() + '(' + tempOutput[0] + ')'
                elif len(tempOutput) == 2:
                    temp = '( ' + tempOutput[0] + ' ' + t.getKid(i).name() + ' ' + tempOutput[1] + ' )'
                else:
                    print('Error. Only binary trees')
                output.append(temp)
                count = counttemp
            else:
                output.append(t.getKid(i).name()) 
                count = count + 1
    if countOld == -1:
        codGPalta = ['plus', 'minus', 'times', 'divide']#, 'cos', 'sin', 'exp', 'sq', 'sqrt']
        codAlgebr = ['+'   , '-'    , '*'    , '/'     ]#, 'cos', 'sin', 'exp', 'sq', 'sqrt']
        for i in range(len(codGPalta)):
            output = output.replace(codGPalta[i], codAlgebr[i])
        return output
    else:
        return output, count