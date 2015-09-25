import pylab
import scipy
import json

class DataReader:
    def __init__(self, fileName):
        self._fileName = fileName

    def read(self, keys):
        data = {}

        with open(self._fileName) as f:

            # title line
            line = f.readline()

            titles = map(lambda s: s.strip(), line.split('\t'))
            keyMap = {}
            for k in keys:
                keyMap[k] = titles.index(k)
                data[k] = []

            # data section
            while True:
                line = f.readline()
                if line == '': break;
                elements = map(lambda s: s.strip(), line.split('\t'))

                for k, idx in keyMap.iteritems():
                    data[k].append(elements[idx])
        return data
    pass

class ConfigReader:
    def __init__(self, fileName):
        self._fileName = fileName
        self._data = None
        pass

    def read(self):
        if self._data == None:
            with open(self._fileName) as data_file:
                self._data = json.load(data_file)
        return self._data

    def getConfigDict(self):
        return self.read()

    def __getitem__(self, k):
        db = self.read()
        return db[k]

    def getDataScale(self, labels):
        db = self.read()['data_scales']
        return [db[i] for i in labels]

class DataPlotter:
    def __init__(self): pass

    def plotNoHVACFromFile(self, fileName, nDays = 14):
        dr = DataReader(fileName)
        data = dr.read(['outdoor', 'space1', 'space2', 'space3', 'space4', 'space5', 'month', 'day', 'hour'])

        self.plotNoHVAC(data, nDays)

    def plotNoHVAC(self, data, nDays = 14):

        pylab.figure(1)
        plt = pylab.subplot(6,1,1)
        plt.vlines(range(24 * nDays), scipy.arange(24 * nDays) * 0, data['outdoor'][0:24*nDays])
        plt.plot(data['outdoor'][0:24*nDays])
        plt.set_xlabel('time (s)')
        plt.set_ylabel('Outdoor (C)')

        for i in scipy.arange(5) + 1:
            plt = pylab.subplot(6, 1, i + 1)
            plt.plot(data['space' + str(i)][0:24*nDays])
            plt.set_xlabel('time (s)')
            plt.set_ylabel('Space #%d (C)' % i)

        pylab.show()

