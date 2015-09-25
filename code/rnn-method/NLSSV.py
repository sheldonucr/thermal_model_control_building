#!/usr/bin/env python

import json
import sys
import os

import numpy as NP
import numpy.random as RND
import pylab as PL

import MyUtils as MU

from scipy import sin, rand, arange

from pybrain.structure                  import FullConnection, RecurrentNetwork, IdentityConnection
from pybrain.structure.modules          import SigmoidLayer, LinearLayer, TanhLayer, BiasUnit
from pybrain.datasets                   import SequentialDataSet
from pybrain.supervised                 import RPropMinusTrainer, BackpropTrainer
from pybrain.tools.validation           import testOnSequenceData, ModuleValidator
from pybrain.tools.shortcuts            import buildNetwork
from pybrain.tools.xml.networkwriter    import NetworkWriter

def main():
    config = MU.ConfigReader('configs/%s' % sys.argv[1])
    config.read()

    logDir = '%s-%s' % (__file__, sys.argv[1])
    os.mkdir(logDir)

    with open('%s/config.txt' % logDir, 'w') as outfile:
        json.dump(config.getConfigDict(), outfile, indent=4)

    dr = MU.DataReader(config['input_tsv_path'])
    data = dr.read(config['interested_columns'])

    inLabels = config['input_columns']

    outLabels = config['output_columns']

    tds, vds = seqDataSetPair(data, inLabels, outLabels, config['seq_label_column'],
            config['test_seqno'], config['validation_seqno'])

    inScale = config.getDataScale(inLabels)
    outScale = config.getDataScale(outLabels)

    normalizeDataSet(tds, ins = inScale, outs = outScale)
    normalizeDataSet(vds, ins = inScale, outs = outScale)

    trainData = tds
    validationData = vds

    fdim = tds.indim / 2 + 15
    xdim = tds.outdim * 2

    rnn = buildNetwork(tds.indim,
            fdim, fdim, xdim,
            tds.outdim,
            hiddenclass=SigmoidLayer,
            recurrent=True)

    rnn.addRecurrentConnection(FullConnection(rnn['hidden2'], rnn['hidden0']))
    rnn.sortModules()

    trainer = RPropMinusTrainer(rnn, dataset=trainData, batchlearning=True, verbose=True)
    #trainer = RPropMinusTrainer(rnn, dataset=trainData, batchlearning=True, verbose=True)
    #trainer = BackpropTrainer(rnn, dataset=trainData, learningrate=0.0001,
    #        lrdecay=1.0, momentum=0.4, verbose=True, batchlearning=False,
    #        weightdecay=0)

    errTime = []
    errTrain = []
    errValidation = []
    epochNo = 0
    while True:

        for i in range(config['epochs_per_update']):
            trainer.train()

        epochNo += config['epochs_per_update']
        NetworkWriter.writeToFile(rnn, '%s/Epoch_%d.xml' % (logDir, epochNo))
        NetworkWriter.writeToFile(rnn, '%s/Latest.xml' % logDir)

        tOut = ModuleValidator.calculateModuleOutput(rnn, trainData)
        vOut = ModuleValidator.calculateModuleOutput(rnn, validationData)

        tScaler = config.getDataScale([config['output_scalar_label']])[0][1]
        tAvgErr = NP.sqrt(NP.mean((trainData['target'] - tOut) ** 2)) * tScaler
        vAvgErr = NP.sqrt(NP.mean((validationData['target'] - vOut) ** 2)) * tScaler

        tMaxErr = NP.max(NP.abs(trainData['target'] - tOut)) * tScaler
        vMaxErr = NP.max(NP.abs(validationData['target'] - vOut)) * tScaler

        errTrain.append(tAvgErr)
        errValidation.append(vAvgErr)
        errTime.append(epochNo)

        print "Training error:      avg %5.3f degC      max %5.3f degC" % (tAvgErr, tMaxErr)
        print "Validation error:    avg %5.3f degC      max %5.3f degC" % (vAvgErr, vMaxErr)
        print "------------------------------------------------------------------------------"

        if (config['visualize_on_training'] == 'yes'):

            PL.figure(1)
            PL.ioff()
            visulizeDataSet(rnn, trainData, 0,
                    config['visualized_columns']['input'],
                    config['visualized_columns']['output'])
            PL.ion()
            PL.draw()

            PL.figure(2)
            PL.ioff()
            visulizeDataSet(rnn, validationData, 0,
                    config['visualized_columns']['input'],
                    config['visualized_columns']['output'])
            PL.ion()
            PL.draw()

            p = PL.figure(3)
            PL.ioff()
            p.clear()
            PL.plot(errTime, errTrain, label = 'Train')
            PL.plot(errTime, errValidation, label = 'Validation')
            PL.legend()
            PL.ion()
            PL.draw()

def addSubNet(nn, prefix, indim, xdim, outdim):

    np = prefix + '_'

    nn.addModule(LinearLayer(indim, name=np+'in'))
    nn.addModule(LinearLayer(outdim, name=np+'out'))

    nn.addModule(SigmoidLayer(indim + xdim, name=np+'f0'))
    nn.addModule(SigmoidLayer(indim + xdim, name=np+'f1'))
    nn.addModule(SigmoidLayer(indim + xdim, name=np+'f2'))
    nn.addModule(SigmoidLayer(indim + xdim, name=np+'f3'))
    nn.addModule(SigmoidLayer(xdim, name=np+'x'))

    nn.addConnection(FullConnection(nn[np+'in'], nn[np+'f0'], outSliceTo=indim))

    nn.addConnection(FullConnection(nn[np+'f0'], nn[np+'f1'], name=np+'f0~f1'))
    nn.addConnection(FullConnection(nn[np+'f1'], nn[np+'f2'], name=np+'f1~f2'))
    nn.addConnection(FullConnection(nn[np+'f2'], nn[np+'f3'], name=np+'f2~f3'))
    nn.addConnection(FullConnection(nn[np+'f3'], nn[np+'x'], name=np+'f3~x'))

    nn.addRecurrentConnection(FullConnection(nn[np+'x'], nn[np+'f0'], outSliceFrom=indim))

    for i in range(outdim):
        nn.addConnection(FullConnection(nn[np+'x'], nn[np+'out'], inSliceFrom=i, inSliceTo=i+1, outSliceFrom=i, outSliceTo=i+1))
    nn.addConnection(FullConnection(nn['b'], nn[np+'out']))

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

def visulizeDataSet(network, data, seqno, in_labels, out_labels):

    seq = data.getSequence(seqno)
    tmpDs = SequentialDataSet(data.indim, data.outdim)
    tmpDs.newSequence()

    for i in xrange(data.getSequenceLength(seqno)):
        tmpDs.addSample(seq[0][i], seq[1][i])

    nplots = len(in_labels) + len(out_labels)

    for i in range(len(in_labels)):
        p = PL.subplot(nplots, 1, i + 1)
        p.clear()
        p.plot(tmpDs['input'][:, i])
        p.set_ylabel(in_labels[i])

    for i in range(len(out_labels)):
        p = PL.subplot(nplots, 1, i + 1 + len(in_labels))
        p.clear()

        output = ModuleValidator.calculateModuleOutput(network, tmpDs)

        p.plot(tmpDs['target'][:, i], label='train')
        p.plot(output[:, i], label='sim')

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
