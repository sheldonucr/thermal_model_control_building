#!/usr/bin/env python

import json
import sys
import os
import re

import numpy as NP
import numpy.random as RND
import pylab as PL

import MyUtils as MU

from scipy import sin, rand, arange

from pybrain.structure.modules          import LSTMLayer, SoftmaxLayer, SigmoidLayer
from pybrain.datasets                   import SequentialDataSet
from pybrain.supervised                 import RPropMinusTrainer, BackpropTrainer
from pybrain.tools.validation           import testOnSequenceData, ModuleValidator
from pybrain.tools.shortcuts            import buildNetwork
from pybrain.tools.xml.networkwriter    import NetworkWriter
from pybrain.tools.xml.networkreader    import NetworkReader

def main():
    logDir = sys.argv[1]
    config = MU.ConfigReader('%s/%s' % (logDir, 'config.txt'))
    config.read()

    dr = MU.DataReader(config['input_tsv_path'])
    data = dr.read(config['interested_columns'])

    inLabels = config['input_columns']

    outLabels = config['output_columns']

    tds, vds = seqDataSetPair(data, inLabels, outLabels,
            config['seq_label_column'], config['test_seqno'],
            config['validation_seqno'])

    inScale = config.getDataScale(inLabels)
    outScale = config.getDataScale(outLabels)

    normalizeDataSet(tds, ins = inScale, outs = outScale)
    normalizeDataSet(vds, ins = inScale, outs = outScale)

    trainData = tds
    validationData = vds

    rnn = NetworkReader.readFrom('%s/%s' % (logDir, 'Latest.xml'))

    tOut = ModuleValidator.calculateModuleOutput(rnn, trainData)
    vOut = ModuleValidator.calculateModuleOutput(rnn, validationData)

    tScaler = config.getDataScale([config['output_scalar_label']])[0][1]
    tAvgErr = NP.sqrt(NP.mean((trainData['target'] - tOut) ** 2)) * tScaler
    vAvgErr = NP.sqrt(NP.mean((validationData['target'] - vOut) ** 2)) * tScaler

    tMaxErr = NP.max(NP.abs(trainData['target'] - tOut)) * tScaler
    vMaxErr = NP.max(NP.abs(validationData['target'] - vOut)) * tScaler

    print "Training error:      avg %5.3f degC      max %5.3f degC" % (tAvgErr, tMaxErr)
    print "Validation error:    avg %5.3f degC      max %5.3f degC" % (vAvgErr, vMaxErr)

    if len(sys.argv) == 3:
        filename = sys.argv[2]
        outbuf = []

        column_titles = [config['seq_label_column'], "training", "validation"]
        column_titles += config['input_columns']
        column_titles += map(lambda s: "%s_actual" % s, config['output_columns'])
        column_titles += map(lambda s: "%s_simulated" % s, config['output_columns'])

        outbuf.append(reduce(lambda lhs, rhs: "%s\t%s" % (lhs, rhs), column_titles))

        for data in [tds, vds]:
            print "Generating %s set curves..." % ("training" if data==tds else "validation")
            for seqno in xrange(data.getNumSequences()):
                print "    %s sequence #%d..." % (("Training" if data==tds else "Validation"), seqno)
                seq = data.getSequence(seqno)
                tmpDs = SequentialDataSet(data.indim, data.outdim)
                tmpDs.newSequence()
                for i in xrange(data.getSequenceLength(seqno)):
                    tmpDs.addSample(seq[0][i], seq[1][i])

                output = ModuleValidator.calculateModuleOutput(rnn, tmpDs)

                # denormalize
                for i in range(len(config['input_columns'])):
                    tmpDs['input'][:, i] *= inScale[i][1]
                    tmpDs['input'][:, i] += inScale[i][0]
                for i in range(len(config['output_columns'])):
                    tmpDs['target'][:, i] *= outScale[i][1]
                    tmpDs['target'][:, i] += outScale[i][0]

                    output[:, i] *= outScale[i][1]
                    output[:, i] += outScale[i][0]

                for i in xrange(data.getSequenceLength(seqno)):

                    line = []
                    line += [seqno, 1 if data==tds else 0, 1 if data==vds else 0]
                    line += tmpDs.getSample(i)[0].tolist()
                    line += tmpDs.getSample(i)[1].tolist()
                    line += output[i].tolist()
                    outbuf.append(reduce(lambda lhs, rhs: "%s\t%s" % (lhs, rhs), line))

                pass # for
            pass # for

        print "Writing results into file '%s'..." % filename
        with open(filename, "w") as f:
            for line in outbuf:
                print >> f, line
        return

    while True:
        s = raw_input("""
Please type what you want to do:
    <N> ------------ Plot the N-th normalized curves in input data set
    <{T|V} N> ------ Plot the N-th normalized curves in Train|Validation set
    <W filename> --- Write the output to file 'filename'
""")
        s = s.strip()
        if s == '': continue

        m = re.match('(\d+)', s)
        if m:
            seqno = int(m.group(1))
            mIndex = None
            mData = None
            mString = None
            if seqno in config['test_seqno']:
                mIndex = config['test_seqno'].index(seqno)
                mData = trainData
                mString = "training"
            elif seqno in config['validation_seqno']:
                mIndex = config['validation_seqno'].index(seqno)
                mData = validationData
                mString = "validation"
            else:
                print "Sequence #%d is not available." % seqno
                continue

            print "Visualizing Sequence #%d (%s data)..." % (seqno, mString)

            PL.figure(1)
            visulizeDataSet(rnn, mData, mIndex,
                    config['visualized_columns']['input'],
                    outLabels, inLabels, outLabels)
            PL.suptitle("Sequence #%d (%s)" % (seqno, mString))
            PL.show()
            continue

        m = re.match('([TVtv])(\d+)', s)
        if m:

            dsString = None
            dsData = None
            dsSeqNO = int(m.group(2))
            dsRealSeq = None

            if m.group(1).upper() == 'T':
                dsString = 'vraining'
                dsData = trainData
                dsRealSeq = config['test_seqno'][dsSeqNO]

            else:
                dsString = 'validation'
                dsData = validationData
                dsRealSeq = config['test_seqno'][dsSeqNO]

            if not dsSeqNO in range(dsData.getNumSequences()):
                print "Sequence #%d is not available." % dsSeqNO
                continue

            print "Visualizing %s sequence #%d (realseq #%d)..." % (dsString, dsSeqNO, dsRealSeq)
            PL.figure(1)
            visulizeDataSet(rnn, dsData, dsSeqNO,
                    config['visualized_columns']['input'], outLabels, inLabels,
                    outLabels)
            PL.suptitle("Sequence #%d in %s dataset (realseq #%d)" % (dsSeqNO, dsString, dsRealSeq))
            visulizeAll5Inputs(rnn, dsData, dsSeqNO, inLabels)
            PL.show()

            continue

        m = re.match('[Ww]\s+(\S+)', s)
        if m:
            filename = "%s.tsv" % m.group(1)
            outbuf = []

            column_titles = [config['seq_label_column'], "training", "validation"]
            column_titles += config['input_columns']
            column_titles += map(lambda s: "%s_actual" % s, config['output_columns'])
            column_titles += map(lambda s: "%s_simulated" % s, config['output_columns'])

            outbuf.append(reduce(lambda lhs, rhs: "%s\t%s" % (lhs, rhs), column_titles))

            for data in [tds, vds]:
                print "Generating %s set curves..." % ("training" if data==tds else "validation")
                for seqno in xrange(data.getNumSequences()):
                    print "    %s sequence #%d..." % (("Training" if data==tds else "Validation"), seqno)
                    seq = data.getSequence(seqno)
                    tmpDs = SequentialDataSet(data.indim, data.outdim)
                    tmpDs.newSequence()
                    for i in xrange(data.getSequenceLength(seqno)):
                        tmpDs.addSample(seq[0][i], seq[1][i])

                    output = ModuleValidator.calculateModuleOutput(rnn, tmpDs)

                    # denormalize
                    for i in range(len(config['input_columns'])):
                        tmpDs['input'][:, i] *= inScale[i][1]
                        tmpDs['input'][:, i] += inScale[i][0]
                    for i in range(len(config['output_columns'])):
                        tmpDs['target'][:, i] *= outScale[i][1]
                        tmpDs['target'][:, i] += outScale[i][0]

                        output[:, i] *= outScale[i][1]
                        output[:, i] += outScale[i][0]

                    for i in xrange(data.getSequenceLength(seqno)):

                        line = []
                        line += [seqno, 1 if data==tds else 0, 1 if data==vds else 0]
                        line += tmpDs.getSample(i)[0].tolist()
                        line += tmpDs.getSample(i)[1].tolist()
                        line += output[i].tolist()
                        outbuf.append(reduce(lambda lhs, rhs: "%s\t%s" % (lhs, rhs), line))

                    pass # for
                pass # for

            print "Writing results into file '%s'..." % filename
            with open(filename, "w") as f:
                for line in outbuf:
                    print >> f, line
            continue # file exporting

        print "Invalid input. Examples: 3, 9, T10, V8, t11, v1, w output"

def normalizeDataSet(data, ins = None, outs = None):
    inscale = []
    outscale = []

    for i in range(data.indim):
        if ins == None:
            mu = NP.mean(data['input'][:, i])
            sigma = NP.std(data['input'][:, i])
        else:
            mu, sigma = ins[i]

        data['input'][:, i] -= mu
        data['input'][:, i] /= sigma

        inscale.append((mu, sigma))

    for i in range(data.outdim):

        if outs == None:

            maxPossible = NP.max(data['target'][:, i])
            minPossible = NP.min(data['target'][:, i])
            mu = minPossible
            sigma = maxPossible - minPossible
        else:
            mu, sigma = outs[i]

        data['target'][:, i] -= mu
        data['target'][:, i] /= sigma

        outscale.append((mu, sigma))

    return (inscale, outscale)

def visulizeAll5Inputs(network, data, seqno, labels):
    figNO = 100
    for i in [1,2,3,4,5]:
        PL.figure(figNO)
        visulizeInput(network, data, seqno, labels, "_%d" % i)
        figNO += 1

def visulizeInput(network, data, seqno, labels, keyword):
    inLabels = []
    for l in labels:
        if re.search(keyword, l):
            inLabels.append(l)

    visulizeDataSet(network, data, seqno, inLabels, [], labels, [])

def visulizeDataSet(network, data, seqno, in_labels, out_labels, in_pool, out_pool):

    seq = data.getSequence(seqno)
    tmpDs = SequentialDataSet(data.indim, data.outdim)
    tmpDs.newSequence()

    for i in xrange(data.getSequenceLength(seqno)):
        tmpDs.addSample(seq[0][i], seq[1][i])

    nplots = len(in_labels) + len(out_labels)

    for i in range(len(in_labels)):
        p = PL.subplot(nplots, 1, i + 1)
        p.clear()
        p.plot(tmpDs['input'][:, in_pool.index(in_labels[i])])
        p.set_ylabel(in_labels[i])

    for i in range(len(out_labels)):
        p = PL.subplot(nplots, 1, i + 1 + len(in_labels))
        p.clear()

        output = ModuleValidator.calculateModuleOutput(network, tmpDs)

        p.plot(tmpDs['target'][:, out_pool.index(out_labels[i])], label='train')
        p.plot(output[:, out_pool.index(out_labels[i])], label='sim')

        p.legend()
        p.set_ylabel(out_labels[i])

def seqDataSetPair(data, in_labels, out_labels, seq_title, tseqs, vseqs):

    tds = SequentialDataSet(len(in_labels), len(out_labels))
    vds = SequentialDataSet(len(in_labels), len(out_labels))
    ds = None

    for i in xrange(len(data[in_labels[0]])):

        if i == 0 or data[seq_title][i] != data[seq_title][i - 1]:
            if int(data[seq_title][i]) in tseqs:
                ds = tds
                ds.newSequence()
            elif int(data[seq_title][i]) in vseqs:
                ds = vds
                ds.newSequence()
            else:
                ds = None

        if ds == None: continue

        din = [data[l][i] for l in in_labels]
        dout = [data[l][i] for l in out_labels]

        ds.addSample(din, dout)

    return (tds, vds)

if __name__ == '__main__':
    main()
